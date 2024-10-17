from typing import Optional, Tuple

from numpy.typing import ArrayLike

from matplotlib import pyplot as plt
from matplotlib.pyplot import Figure, Axes
from sklearn.metrics import r2_score
from scipy.stats import spearmanr


def plot_prediction_scatter(
    y_true: ArrayLike,
    y_pred: ArrayLike,
    x_label: str = "Ground truth",
    y_label: str = "Prediction",
    title: Optional[str] = None,
    axis_label_fontsize: int = 16,
    title_fontsize: int = 18,
    margin_ratio: float = 0.15,
) -> Tuple[Figure, Axes]:
    x_min, x_max = min(y_true), max(y_true)
    y_min, y_max = min(y_pred), max(y_pred)
    x_margin = margin_ratio * (x_max - x_min)
    y_margin = margin_ratio * (y_max - y_min)
    fig, ax = plt.subplots()
    ax.scatter(y_true, y_pred, c="k", s=50, alpha=0.75)
    ax.set_xlabel(x_label, fontsize=axis_label_fontsize)
    ax.set_ylabel(y_label, fontsize=axis_label_fontsize)
    ax.set_title(title, fontsize=title_fontsize)
    ax.set_xlim(x_min, x_max + x_margin)
    ax.set_ylim(y_min, y_max + y_margin)
    spearman, _ = spearmanr(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    metrics_str = f"Spearman: {spearman:.2f}\nR2: {r2:.2f}"
    box_format = dict(boxstyle="round", facecolor="white", alpha=0.5)
    ax.text(
        0.05,
        0.95,
        metrics_str,
        transform=ax.transAxes,
        fontsize=12,
        verticalalignment="top",
        bbox=box_format,
    )
    return fig, ax
