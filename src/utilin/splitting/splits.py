import pandas as pd
import numpy as np


def assign_random_folds(
    data: pd.DataFrame, n_folds: int = 5, seed: int = None
) -> pd.DataFrame:
    if seed is not None:
        np.random.seed(seed)
    df = data.copy()
    fold_col = f"fold_random_{n_folds}"
    folds = np.tile(np.arange(n_folds), np.ceil(len(df) / n_folds).astype(int))[
        : len(df)
    ]
    np.random.shuffle(folds)
    df[fold_col] = folds
    return df


def assign_modulo_folds(
    data: pd.DataFrame,
    mutation_col: str = "mutant",
    n_folds: int = 5,
    fold_start: int = 0,
) -> pd.DataFrame:
    df = data.copy()
    fold_col = f"fold_modulo_{n_folds}"
    mutated_positions = df[mutation_col].map(lambda x: int(x[1:-1])).values
    folds = (mutated_positions - 1 + fold_start) % n_folds
    df[fold_col] = folds
    return df


def assign_contiguous_folds(
    data: pd.DataFrame,
    mutation_col: str = "mutant",
    n_folds: int = 5,
    span: int = 5,
    fold_start: int = 0,
) -> pd.DataFrame:
    sequence_length = len(data["sequence"][0])
    assert (
        np.ceil(sequence_length / n_folds) >= span
    ), f"Span is too large for a sequence of length {sequence_length} and {n_folds} folds"
    df = data.copy()
    fold_col = f"fold_contiguous_{n_folds}"
    n_regions_per_fold = np.ceil(sequence_length / (span * n_folds)).astype(int)
    mutated_positions = df[mutation_col].map(lambda x: int(x[1:-1])).values
    position_fold_ids = np.tile(
        np.repeat(np.arange(n_folds) + fold_start, span) % n_folds, n_regions_per_fold
    )[:sequence_length]
    folds = list(map(lambda x: position_fold_ids[x - 1], mutated_positions))
    df[fold_col] = folds
    return df
