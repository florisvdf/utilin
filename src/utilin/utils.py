import requests
from pathlib import Path
from typing import Union

from biotite.sequence.io.fasta import FastaFile


def fetch_uniprot_sequence(uniprot_id: str) -> str:
    response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
    body = response.text
    sequence = "".join(body.split("\n")[1:])
    return sequence


def read_fasta(path: Union[str, Path]) -> FastaFile:
    return FastaFile.read(str(path))


def sequence_to_mutations(variant: str, reference: str) -> str:
    return ":".join(
        [
            f"{aa_ref}{pos+1}{aa_var}"
            for pos, (aa_ref, aa_var) in enumerate(zip(reference, variant))
            if aa_ref != aa_var
        ]
    )
