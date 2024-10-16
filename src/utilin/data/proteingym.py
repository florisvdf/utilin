import re
from pathlib import Path

import pandas as pd
import requests
from biotite.sequence import ProteinSequence
from biotite.sequence.align import SubstitutionMatrix, align_optimal
from loguru import logger
from polyleven import levenshtein

from utilin.utils import fetch_uniprot_sequence, read_fasta


class ProteinGym:
    """
    Currently only support for substitutions.
    """

    def __init__(self, proteingym_location, meta_data_path):
        self.proteingym_location = proteingym_location
        self.meta_data_path = meta_data_path
        self.reference_information = pd.read_csv(self.meta_data_path, index_col=None)
        self.available_pdbs = {}

    def update_reference_information(self):
        logger.info("Updating reference information with structure information.")
        pdb_entry_pattern = r".*PDB; [A-Z0-9]{4};.*"
        region_pattern = r"A=(\d+-\d+)"
        for i, row in self.reference_information.iterrows():
            uniprot_id = row["UniProt_ID"]
            region_mutated = row["region_mutated"]
            response = requests.get(
                f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.txt"
            )
            if response.status_code == 400:
                continue
            body = response.text
            pdb_entries = list(
                filter(
                    lambda x: True if re.match(pdb_entry_pattern, x) else False,
                    body.split("\n"),
                )
            )
            reference_sequence = self.reference_information.loc[i, "target_seq"]
            uniprot_sequence = fetch_uniprot_sequence(uniprot_id)
            self.available_pdbs[uniprot_id] = pdb_entries
            self.reference_information.loc[i, "has_pdb_structure"] = (
                len(pdb_entries) > 0
            )
            self.reference_information.loc[i, "structure_covers_mutated_region"] = False
            self.reference_information.loc[i, "uniprot_sequence"] = uniprot_sequence
            self.reference_information.loc[i, "reference_distance_to_uniprot"] = (
                self.distance_of_reference_to_uniprot(
                    reference_sequence, uniprot_sequence
                )
            )
            for entry in pdb_entries:
                match = re.search(region_pattern, entry)
                if match:
                    structure_region = match.group(1)
                    structure_covers_mutated_region = self.region_is_subregion(
                        region_mutated, structure_region
                    )
                    if structure_covers_mutated_region:
                        self.reference_information.loc[
                            i, "structure_covers_mutated_region"
                        ] = True
                        break

    def describe_dataset(self, dataset_name):
        return self.reference_information[
            self.reference_information["DMS_id"] == dataset_name
        ]

    def fetch_msa(self, dataset_name):
        uniprot_id = self.describe_dataset(dataset_name)["UniProt_ID"].values[0]
        matching_msa_paths = list(
            (Path(self.proteingym_location) / "MSA_files/DMS").rglob(f"{uniprot_id}*")
        )
        logger.info(f"Found {len(matching_msa_paths)} matching MSA files")
        return read_fasta(matching_msa_paths[0])

    def prepare_dataset(self, dataset_name):
        data = pd.read_csv(
            Path(self.proteingym_location)
            / f"cv_folds_singles_substitutions/{dataset_name}.csv"
        ).rename(columns={"mutated_sequence": "sequence"})
        return data

    def distance_of_reference_to_uniprot(self, reference: str, uniprot: str):
        mat = SubstitutionMatrix.std_protein_matrix()
        reference_sequence = ProteinSequence(reference)
        uniprot_sequence = ProteinSequence(uniprot)
        alignment = align_optimal(reference_sequence, uniprot_sequence, mat)[0]
        trace = alignment.trace
        mask = trace[:, 0] != -1
        relevant_uniprot_region = trace[mask][:, 1]
        relevant_uniprot_sequence = uniprot_sequence[relevant_uniprot_region]
        edit_distance = levenshtein(reference, relevant_uniprot_sequence.__str__())
        return edit_distance

    def region_is_subregion(self, region1: str, region2: str):
        start1, end1 = tuple(int(value) for value in region1.split("-"))
        start2, end2 = tuple(int(value) for value in region2.split("-"))
        return (start1 >= start2) & (end1 <= end2)
