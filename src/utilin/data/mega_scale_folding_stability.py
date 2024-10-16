import re

import pandas as pd


class MegaScaleProteinFoldingStability:
    def __init__(self, dataset_path: str):
        self.dataset = pd.read_csv(dataset_path)
        self.dataset["n_mutations"] = self.dataset["mut_type"].map(
            self.find_number_of_mutations
        )
        self.dataset["mutated_positions"] = self.dataset["mut_type"].map(
            self.find_mutated_positions
        )
        self.dataset["synonymous"] = self.dataset["mut_type"].map(self.is_synonymous)
        self.dataset["is_ins_or_del"] = self.dataset["mut_type"].map(
            self.is_insertion_or_deletion
        )
        self.dataset["sequence"] = self.dataset["aa_seq"]

    @property
    def dataset_names(self):
        return self.dataset["WT_name"].map(lambda x: x.split(".")[0]).unique()

    @staticmethod
    def find_number_of_mutations(mut_type):
        if mut_type in ["wt", ""]:
            return 0
        if re.search(r"ins|del", mut_type):
            return 1
        else:
            return len(mut_type.split(":"))

    @staticmethod
    def find_mutated_positions(mut_type):
        return (
            frozenset(int(mut[1:-1]) for mut in mut_type.split(":"))
            if not re.search(r"wt|ins|del", mut_type)
            else frozenset()
        )

    @staticmethod
    def is_synonymous(mut_type):
        return (
            all([mut[0] == mut[-1] for mut in mut_type.split(":")])
            if not re.search(r"wt|ins|del", mut_type)
            else False
        )

    @staticmethod
    def is_insertion_or_deletion(mut_type):
        return True if re.search(r"ins|del", mut_type) else False

    def fetch_dataset(
        self,
        pdb_id,
        keep_synonymous=False,
        keep_insertions_and_deletions=False,
        drop_duplicate_sequences=True,
        order="second",
    ):
        return self.dataset[
            (self.dataset["WT_name"] == f"{pdb_id}.pdb")
            & (self.dataset["synonymous"] == (False or keep_synonymous))
            & (
                self.dataset["is_ins_or_del"]
                == (False or keep_insertions_and_deletions)
            )
            & (
                self.dataset["n_mutations"].isin(
                    {"first": [1], "second": [2], "both": [1, 2]}[order]
                )
            )
        ].drop_duplicates(subset="sequence" if drop_duplicate_sequences else None)

    def get_reference_sequence(self, pdb_id):
        return self.dataset[self.dataset["name"] == f"{pdb_id}.pdb"]["sequence"]

    def describe_dataset(self, pdb_id):
        raise NotImplementedError
