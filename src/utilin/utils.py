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


def mutations_to_sequence(mutations: List[str], reference: str) -> str:
    split_sequence = list(reference)
    for mutation in mutations:
        position = int(mutation[1:-1])
        mutated_residue = mutation[-1]
        split_sequence[position - 1] = mutated_residue
    return "".join(split_sequence)


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


def sets_share_positions(set_1: tuple, set_2: tuple) -> bool:
    return any([x == y for x in set_1 for y in set_2])


def group_position_non_overlapping(combinations: List[tuple]) -> Dict[tuple, int]:
    initial_grouping = dict(Counter(combinations))
    combinations_have_overlap = True
    new_grouping = initial_grouping.copy()
    while combinations_have_overlap:
        unique_sets = list(new_grouping.keys())
        overlap_indicator = np.zeros(
            (len(unique_sets), len(unique_sets)), dtype=np.int32
        )
        for i, set_1 in enumerate(unique_sets):
            for j, set_2 in enumerate(unique_sets):
                overlap_indicator[i, j] = int(sets_share_positions(set_1, set_2))
        combinations_have_overlap = np.sum(np.triu(overlap_indicator, k=1)) > 0
        set_with_most_overlap_idx = np.argmax(np.sum(overlap_indicator, axis=0))
        sets_to_combine = [
            unique_sets[r]
            for r in np.where(overlap_indicator[set_with_most_overlap_idx])[0]
        ]
        new_combination = tuple(
            set(position for combination in sets_to_combine for position in combination)
        )
        counts = 0
        for old_set in sets_to_combine:
            counts += new_grouping[old_set]
            new_grouping.pop(old_set, None)
        new_grouping[new_combination] = counts
    return new_grouping
