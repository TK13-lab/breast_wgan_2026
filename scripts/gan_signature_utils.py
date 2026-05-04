from __future__ import annotations

from dataclasses import dataclass
import os
from pathlib import Path
import random
import re
import warnings

MPL_CACHE_DIR = Path(__file__).resolve().parents[1] / "cache" / "matplotlib"
MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(MPL_CACHE_DIR))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, average_precision_score, f1_score, roc_auc_score
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
import torch
import torch.nn as nn
from torch.utils.data import DataLoader, TensorDataset


LABEL_PATTERN = re.compile(r"^(BN1|BN2|DN|EN|DT|ET|FT|GT|BT)")
LABEL_MAP = {
    "BN1": 0,
    "BN2": 0,
    "DN": 0,
    "EN": 0,
    "DT": 1,
    "ET": 1,
    "FT": 1,
    "GT": 1,
    "BT": 1,
}


@dataclass
class TrainConfig:
    n_epochs: int = 50
    batch_size: int = 128
    lr: float = 1e-3
    weight_decay: float = 1e-5
    alpha: float = 1.0
    beta: float = 1.0
    gamma: float = 0.1
    lambda_gp: float = 10.0
    n_critic: int = 3
    device: str = "cpu"


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)


def resolve_device(requested: str = "auto") -> str:
    requested = requested.lower()
    if requested == "auto":
        return "cuda" if torch.cuda.is_available() else "cpu"
    if requested == "cuda" and not torch.cuda.is_available():
        warnings.warn("CUDA requested but not available. Falling back to CPU.", RuntimeWarning)
        return "cpu"
    return requested


def manual_stratified_split(df_t: pd.DataFrame, seed: int, test_size: float) -> tuple[pd.DataFrame, pd.DataFrame]:
    train_chunks = []
    test_chunks = []

    for label, group in df_t.groupby("label", sort=True, group_keys=False):
        n_obs = len(group)
        if n_obs < 2:
            raise ValueError(f"Class {label} has only {n_obs} sample(s); a holdout split is not possible.")

        lower = 1.0 / n_obs + 1e-9
        upper = 1.0 - 1.0 / n_obs - 1e-9
        test_size_eff = min(max(test_size, lower), upper)

        n_test = max(1, int(round(n_obs * test_size_eff)))
        if n_test >= n_obs:
            n_test = n_obs - 1

        shuffled = group.sample(frac=1.0, random_state=seed).reset_index(drop=True)
        test_chunks.append(shuffled.iloc[:n_test].copy())
        train_chunks.append(shuffled.iloc[n_test:].copy())

    train_df = pd.concat(train_chunks, axis=0).sample(frac=1.0, random_state=seed).reset_index(drop=True)
    test_df = pd.concat(test_chunks, axis=0).sample(frac=1.0, random_state=seed).reset_index(drop=True)
    return train_df, test_df


def read_and_prepare_inputs(
    raw_counts_csv: Path,
    seed: int,
    test_size: float,
    exclude_prefixes: tuple[str, ...] = ("RPS", "RPL"),
) -> dict[str, object]:
    df = pd.read_csv(raw_counts_csv, index_col=0, low_memory=False)
    df_t = df.T.reset_index().rename(columns={"index": "sample_id"})

    df_t["label"] = df_t["sample_id"].str.extract(LABEL_PATTERN.pattern, expand=False).map(LABEL_MAP)
    df_t = df_t.dropna(subset=["label"]).copy()
    df_t["label"] = df_t["label"].astype(int)

    train_df, test_df = manual_stratified_split(df_t, seed=seed, test_size=test_size)

    train_features = train_df.drop(columns=["sample_id"]).copy()
    test_features = test_df.drop(columns=["sample_id"]).copy()

    train_gene_cols = [
        col for col in train_features.columns
        if col != "label" and not any(col.startswith(prefix) for prefix in exclude_prefixes)
    ]
    test_gene_cols = [
        col for col in test_features.columns
        if col != "label" and not any(col.startswith(prefix) for prefix in exclude_prefixes)
    ]
    common_gene_cols = sorted(set(train_gene_cols) & set(test_gene_cols))
    if not common_gene_cols:
        raise ValueError(f"No common gene columns remained after filtering for {raw_counts_csv.name}.")

    x_train = train_features[common_gene_cols].copy()
    x_test = test_features[common_gene_cols].copy()
    y_train = train_features["label"].astype(int).to_numpy()
    y_test = test_features["label"].astype(int).to_numpy()

    scaler = StandardScaler(with_mean=True, with_std=True)
    x_train_scaled = scaler.fit_transform(x_train.values.astype(np.float32))
    x_test_scaled = scaler.transform(x_test.values.astype(np.float32))

    x_train_df = pd.DataFrame(x_train_scaled, columns=common_gene_cols)
    x_test_df = pd.DataFrame(x_test_scaled, columns=common_gene_cols)

    return {
        "x_train_df": x_train_df,
        "x_test_df": x_test_df,
        "y_train": y_train,
        "y_test": y_test,
        "train_sample_ids": train_df["sample_id"].tolist(),
        "test_sample_ids": test_df["sample_id"].tolist(),
        "metadata": {
            "n_cells_total": int(df_t.shape[0]),
            "n_genes_raw": int(df.shape[0]),
            "n_genes_model": int(len(common_gene_cols)),
            "n_train": int(x_train_df.shape[0]),
            "n_test": int(x_test_df.shape[0]),
            "train_label_0": int((y_train == 0).sum()),
            "train_label_1": int((y_train == 1).sum()),
            "test_label_0": int((y_test == 0).sum()),
            "test_label_1": int((y_test == 1).sum()),
        },
    }


def propose_architecture(n_input: int) -> tuple[list[int], int]:
    if n_input <= 500:
        return [256, 128], 64
    if n_input <= 3000:
        return [512, 256], 128
    return [1024, 512], 256


def mlp_stack(dims: list[int], activation: str = "gelu", dropout: float = 0.2, layernorm: bool = True) -> nn.Sequential:
    layers: list[nn.Module] = []
    activation_layer: nn.Module = nn.GELU() if activation == "gelu" else nn.ReLU()

    for idx in range(len(dims) - 1):
        layers.append(nn.Linear(dims[idx], dims[idx + 1]))
        if layernorm:
            layers.append(nn.LayerNorm(dims[idx + 1]))
        # This mirrors the legacy notebook, which also applies activation and
        # optional dropout after the final linear layer.
        layers.append(activation_layer.__class__())
        if dropout and dropout > 0:
            layers.append(nn.Dropout(dropout))

    return nn.Sequential(*layers)


class AAEWGANGP(nn.Module):
    def __init__(self, input_dim: int, enc_hidden: list[int], latent_dim: int, cls_hidden: list[int]) -> None:
        super().__init__()
        self.encoder = mlp_stack([input_dim] + enc_hidden + [latent_dim], activation="gelu", dropout=0.2, layernorm=True)
        self.decoder = mlp_stack([latent_dim] + enc_hidden[::-1] + [input_dim], activation="gelu", dropout=0.0, layernorm=False)
        self.classifier = mlp_stack([latent_dim] + cls_hidden + [2], activation="gelu", dropout=0.2, layernorm=True)
        self.disc = nn.Sequential(
            nn.utils.spectral_norm(nn.Linear(latent_dim, 128)),
            nn.GELU(),
            nn.utils.spectral_norm(nn.Linear(128, 32)),
            nn.GELU(),
            nn.utils.spectral_norm(nn.Linear(32, 1)),
        )

    def encode(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)

    def decode(self, z: torch.Tensor) -> torch.Tensor:
        return self.decoder(z)

    def classify(self, z: torch.Tensor) -> torch.Tensor:
        return self.classifier(z)

    def D(self, z: torch.Tensor) -> torch.Tensor:
        return self.disc(z)

    def forward(self, x: torch.Tensor) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        z = self.encode(x)
        x_hat = self.decode(z)
        logits = self.classify(z)
        d_score = self.D(z)
        return z, x_hat, logits, d_score


def make_loaders(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_valid: np.ndarray,
    y_valid: np.ndarray,
    batch_size: int,
) -> tuple[DataLoader, DataLoader]:
    train_ds = TensorDataset(torch.tensor(x_train, dtype=torch.float32), torch.tensor(y_train, dtype=torch.long))
    valid_ds = TensorDataset(torch.tensor(x_valid, dtype=torch.float32), torch.tensor(y_valid, dtype=torch.long))
    return (
        DataLoader(train_ds, batch_size=batch_size, shuffle=True),
        DataLoader(valid_ds, batch_size=batch_size, shuffle=False),
    )


def grad_penalty(discriminator: nn.Module, real_z: torch.Tensor, fake_z: torch.Tensor, cfg: TrainConfig) -> torch.Tensor:
    eps = torch.rand(real_z.size(0), 1, device=cfg.device)
    eps = eps.expand_as(real_z)
    interp = eps * real_z + (1 - eps) * fake_z
    interp.requires_grad_(True)
    d_interp = discriminator(interp)
    grads = torch.autograd.grad(
        outputs=d_interp,
        inputs=interp,
        grad_outputs=torch.ones_like(d_interp),
        create_graph=True,
        retain_graph=True,
        only_inputs=True,
    )[0]
    return ((grads.norm(2, dim=1) - 1.0) ** 2).mean()


def _safe_roc_auc(y_true: np.ndarray, y_score: np.ndarray) -> float:
    return float(roc_auc_score(y_true, y_score)) if len(np.unique(y_true)) > 1 else float("nan")


def _safe_average_precision(y_true: np.ndarray, y_score: np.ndarray) -> float:
    return float(average_precision_score(y_true, y_score)) if len(np.unique(y_true)) > 1 else float("nan")


@torch.no_grad()
def evaluate(model: AAEWGANGP, loader: DataLoader, cfg: TrainConfig) -> tuple[float, float, float, float]:
    model.eval()
    all_probs = []
    all_true = []
    total_loss = 0.0
    total_acc = 0.0
    n_obs = 0
    mse = nn.MSELoss()
    ce = nn.CrossEntropyLoss()

    for xb, yb in loader:
        xb = xb.to(cfg.device)
        yb = yb.to(cfg.device)
        z = model.encode(xb)
        x_hat = model.decode(z)
        logits = model.classify(z)

        loss = cfg.alpha * mse(x_hat, xb) + cfg.beta * ce(logits, yb)
        total_loss += loss.item() * xb.size(0)

        probs = torch.softmax(logits, dim=1)[:, 1]
        preds = torch.argmax(logits, dim=1)
        total_acc += float((preds == yb).sum().item())
        n_obs += xb.size(0)
        all_probs.append(probs.detach().cpu().numpy())
        all_true.append(yb.detach().cpu().numpy())

    probs_arr = np.concatenate(all_probs) if all_probs else np.array([])
    true_arr = np.concatenate(all_true) if all_true else np.array([])

    return (
        total_loss / max(n_obs, 1),
        total_acc / max(n_obs, 1),
        _safe_roc_auc(true_arr, probs_arr),
        _safe_average_precision(true_arr, probs_arr),
    )


def train_cv_return_best(
    x: np.ndarray,
    y: np.ndarray,
    enc_hidden: list[int],
    latent_dim: int,
    cls_hidden: list[int],
    cfg: TrainConfig,
    n_splits: int,
    seed: int,
) -> tuple[AAEWGANGP, pd.DataFrame]:
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    best_model: AAEWGANGP | None = None
    best_metric = -np.inf
    history_rows: list[dict[str, float | int]] = []

    for fold, (train_idx, valid_idx) in enumerate(skf.split(x, y)):
        x_train, x_valid = x[train_idx], x[valid_idx]
        y_train, y_valid = y[train_idx], y[valid_idx]

        model = AAEWGANGP(
            input_dim=x.shape[1],
            enc_hidden=enc_hidden,
            latent_dim=latent_dim,
            cls_hidden=cls_hidden,
        ).to(cfg.device)

        opt_main = torch.optim.Adam(
            list(model.encoder.parameters()) + list(model.decoder.parameters()) + list(model.classifier.parameters()),
            lr=cfg.lr,
            weight_decay=cfg.weight_decay,
            betas=(0.5, 0.9),
        )
        opt_disc = torch.optim.Adam(
            model.disc.parameters(),
            lr=cfg.lr,
            weight_decay=cfg.weight_decay,
            betas=(0.5, 0.9),
        )

        mse = nn.MSELoss()
        ce = nn.CrossEntropyLoss()
        train_loader, valid_loader = make_loaders(x_train, y_train, x_valid, y_valid, cfg.batch_size)

        best_state = None
        best_fold_metric = -np.inf

        for epoch in range(cfg.n_epochs):
            model.train()
            epoch_loss = 0.0
            n_seen = 0

            for xb, yb in train_loader:
                xb = xb.to(cfg.device)
                yb = yb.to(cfg.device)

                for _ in range(cfg.n_critic):
                    z_fake = model.encode(xb).detach()
                    z_real = torch.randn_like(z_fake)

                    d_fake = model.D(z_fake)
                    d_real = model.D(z_real)
                    gp = grad_penalty(model.D, z_real, z_fake, cfg)
                    loss_disc = -(d_real.mean() - d_fake.mean()) + cfg.lambda_gp * gp

                    opt_disc.zero_grad()
                    loss_disc.backward()
                    opt_disc.step()

                z = model.encode(xb)
                x_hat = model.decode(z)
                logits = model.classify(z)
                gen_adv = -model.D(z).mean()

                loss = cfg.alpha * mse(x_hat, xb) + cfg.beta * ce(logits, yb) + cfg.gamma * gen_adv

                opt_main.zero_grad()
                loss.backward()
                opt_main.step()

                epoch_loss += loss.item() * xb.size(0)
                n_seen += xb.size(0)

            valid_loss, valid_acc, valid_auroc, valid_auprc = evaluate(model, valid_loader, cfg)
            history_rows.append(
                {
                    "fold": fold,
                    "epoch": epoch,
                    "train_loss": epoch_loss / max(n_seen, 1),
                    "val_loss": valid_loss,
                    "val_acc": valid_acc,
                    "val_auroc": valid_auroc,
                    "val_auprc": valid_auprc,
                }
            )

            compare_metric = valid_auroc if np.isfinite(valid_auroc) else -np.inf
            if compare_metric > best_fold_metric:
                best_fold_metric = compare_metric
                best_state = {name: tensor.detach().cpu().clone() for name, tensor in model.state_dict().items()}

        if best_state is None:
            best_state = {name: tensor.detach().cpu().clone() for name, tensor in model.state_dict().items()}

        if best_fold_metric > best_metric:
            best_metric = best_fold_metric
            best_model = AAEWGANGP(
                input_dim=x.shape[1],
                enc_hidden=enc_hidden,
                latent_dim=latent_dim,
                cls_hidden=cls_hidden,
            ).to(cfg.device)
            best_model.load_state_dict(best_state)

    if best_model is None:
        raise RuntimeError("Training did not produce a best model.")

    return best_model, pd.DataFrame(history_rows)


def collect_linear_weights(seq: nn.Sequential) -> list[np.ndarray]:
    return [module.weight.detach().cpu().numpy() for module in seq if isinstance(module, nn.Linear)]


def product_chain(weights: list[np.ndarray]) -> np.ndarray:
    if not weights:
        raise RuntimeError("No linear layers were found for weight chaining.")
    acc = weights[0]
    for idx in range(1, len(weights)):
        acc = weights[idx] @ acc
    return acc


def compute_total_weight(model: AAEWGANGP) -> np.ndarray:
    encoder_weights = collect_linear_weights(model.encoder)
    classifier_weights = collect_linear_weights(model.classifier)
    encoder_chain = product_chain(encoder_weights)
    classifier_chain = product_chain(classifier_weights)
    return classifier_chain @ encoder_chain


def rank_genes(
    weight_matrix: np.ndarray,
    gene_names: list[str],
    disease_index: int = 1,
    use_abs: bool = True,
) -> pd.DataFrame:
    scores = weight_matrix[disease_index, :]
    if use_abs:
        scores = np.abs(scores)
    order = np.argsort(scores)[::-1]
    return pd.DataFrame(
        {
            "gene": np.asarray(gene_names)[order],
            "score": scores[order],
            "rank": np.arange(1, len(order) + 1),
        }
    )


def unique_effective_ks(ks: tuple[int, ...], n_features: int) -> list[int]:
    seen: set[int] = set()
    resolved: list[int] = []
    for k in ks:
        k_eff = min(int(k), int(n_features))
        if k_eff not in seen:
            resolved.append(k_eff)
            seen.add(k_eff)
    return resolved


def evaluate_k_curve(
    rank_df: pd.DataFrame,
    x_train_df: pd.DataFrame,
    y_train: np.ndarray,
    x_test_df: pd.DataFrame,
    y_test: np.ndarray,
    ks: tuple[int, ...],
    clf_kind: str,
    seed: int,
) -> pd.DataFrame:
    ordered_genes = rank_df["gene"].tolist()
    results = []

    for k_eff in unique_effective_ks(ks, x_train_df.shape[1]):
        genes_k = ordered_genes[:k_eff]
        x_train = x_train_df[genes_k].values
        x_test = x_test_df[genes_k].values

        if clf_kind == "knn":
            n_neighbors = min(5, max(1, x_train.shape[0]))
            clf = KNeighborsClassifier(n_neighbors=n_neighbors)
        else:
            clf = LogisticRegression(
                solver="liblinear",
                penalty="l2",
                max_iter=2000,
                class_weight="balanced",
                random_state=seed,
            )

        clf.fit(x_train, y_train)
        probs = clf.predict_proba(x_test)[:, 1]
        preds = (probs >= 0.5).astype(int)

        results.append(
            {
                "k": k_eff,
                "AUROC": _safe_roc_auc(y_test, probs),
                "AUPRC": _safe_average_precision(y_test, probs),
                "F1": float(f1_score(y_test, preds)),
                "ACC": float(accuracy_score(y_test, preds)),
            }
        )

    return pd.DataFrame(results)


def plot_metric_curve(df_curve: pd.DataFrame, metric: str, out_png: Path, title_prefix: str) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(6, 4))
    plt.plot(df_curve["k"], df_curve[metric], marker="o")
    plt.xlabel("Top-k genes")
    plt.ylabel(metric)
    plt.title(f"{title_prefix}: {metric} vs top-k")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
