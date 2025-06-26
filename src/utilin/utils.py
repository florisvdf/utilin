import requests
from pathlib import Path
from collections import defaultdict, Counter
from typing import Union, Dict, List

import numpy as np

from biotite.sequence.io.fasta import FastaFile


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
    body = response.text
    sequence = "".join(body.split("\n")[1:])
    return sequence


def read_fasta(path: Union[str, Path]) -> FastaFile:
    return FastaFile.read(str(path))


def write_fasta(
    path: Union[str, Path], sequences: Union[str, FastaFile, Dict, List]
) -> None:
    file = FastaFile()
    path = str(path)
    if isinstance(sequences, str):
        sequences = read_fasta(sequences)
    if isinstance(sequences, list):
        for i, sequence in enumerate(sequences):
            file[str(i)] = sequence
    elif isinstance(sequences, dict):
        for key, sequence in sequences.items():
            file[key] = sequence
    elif isinstance(sequences, FastaFile):
        file = sequences
    else:
        raise TypeError("Sequences are not of type FastaFile, Dict or List")
    file.write(path)


def reference_from_variant_sequences(sequences: list) -> str:
    lengths = [len(seq) for seq in sequences]
    unique_lengths = np.unique(lengths)
    if len(unique_lengths) != 1:
        raise ValueError("Sequence must all be of the same length!")
    length = lengths[0]
    sequence_grid = np.vstack(list(map(list, sequences)))
    most_frequent = [
        Counter(sequence_grid[:, i]).most_common()[0][0] for i in range(length)
    ]
    reference = "".join(most_frequent)
    return reference


def variant_sequence_to_mutations(variant: str, reference: str) -> str:
    return ":".join(
        [
            f"{aa_ref}{pos+1}{aa_var}"
            for pos, (aa_ref, aa_var) in enumerate(zip(reference, variant))
            if aa_ref != aa_var
        ]
    )


def determine_residue_offset(structure_sequence: str, reference_sequence: str) -> int:
    structure_residue_positions = {
        index: aa for index, aa in enumerate(structure_sequence)
    }
    reference_residue_positions = {
        index: aa for index, aa in enumerate(reference_sequence)
    }
    matching_pairs = []
    for idx_a, char_a in structure_residue_positions.items():
        for idx_b, char_b in reference_residue_positions.items():
            if char_a == char_b:
                matching_pairs.append((idx_a, idx_b))
    offsets = defaultdict(int)
    for idx_a, idx_b in matching_pairs:
        offset = idx_b - idx_a
        offsets[offset] += 1
    if offsets:
        most_frequent_offset = max(offsets, key=offsets.get)
        return most_frequent_offset
    else:
        return 0
