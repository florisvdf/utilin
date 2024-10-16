from pathlib import Path
from typing import Union

import requests
from biotite.sequence.io.fasta import FastaFile


def fetch_uniprot_sequence(uniprot_id: str):
    response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.fasta")
    body = response.text
    sequence = "".join(body.split("\n")[1:])
    return sequence

def read_fasta(path: Union[str, Path]):
    return FastaFile.read(str(path))
