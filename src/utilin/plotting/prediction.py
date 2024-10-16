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
) -> Tuple[Figure, Axes]:
    fig, ax = plt.subplots()
    ax.scatter(y_true, y_pred, c="k", s=50, alpha=0.75)
    ax.set_xlabel(x_label, fontsize=axis_label_fontsize)
    ax.set_ylabel(y_label, fontsize=axis_label_fontsize)
    ax.set_title(title)
    _, spearman = spearmanr(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    metrics_str = f"Spearman: {spearman:.2f}\nR2: {r2:.2f}"
    box_format = dict(boxstyle="round", facecolor="white", alpha=0.5)
    ax.text(
        0.05,
        0.95,
        metrics_str,
        transform=ax.transAxes,
        fontsize=14,
        verticalalignment="top",
        bbox=box_format,
    )
    return fig, ax
