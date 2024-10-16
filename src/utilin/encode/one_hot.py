from typing import Union, List

import pandas as pd
import numpy as np

from utilin.constants import AA_ALPHABET, AA_ALPHABET_GREMLIN


def encode_sequences_one_hot(
    sequences: Union[List[str], pd.Series], aa_alphabet: str = "default"
) -> np.ndarray:
    alphabet = {
        "default": AA_ALPHABET,
        "gremlin": AA_ALPHABET_GREMLIN,
    }[aa_alphabet]
    return np.array(
        [
            np.array([np.eye(len(alphabet))[alphabet.index(aa)] for aa in sequence])
            for sequence in sequences
        ]
    )
