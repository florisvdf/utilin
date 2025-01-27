from typing import Union, List

import pandas as pd
import numpy as np

from utilin.constants import AA_ALPHABET, AA_ALPHABET_GREMLIN, AA_ALPHABET_ESM_IF1


def encode_sequences_one_hot(
    sequences: Union[List[str], pd.Series, str], aa_alphabet: str = "default"
) -> np.ndarray:
    alphabet = {
        "default": AA_ALPHABET,
        "gremlin": AA_ALPHABET_GREMLIN,
        "esm_if1": AA_ALPHABET_ESM_IF1,
    }[aa_alphabet]
    if isinstance(sequences, str):
        sequences = [sequences]
    return np.array(
        [
            np.array([np.eye(len(alphabet))[alphabet.index(aa)] for aa in sequence])
            for sequence in sequences
        ]
    )
