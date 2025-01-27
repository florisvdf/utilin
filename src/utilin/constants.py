AA_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
AA_ALPHABET_GREMLIN = list("ARNDCQEGHILKMFPSTWYV-")
AA_ALPHABET_ESM_IF1 = "LAGVSERTIDPKQNFYMHWC"
THREE_TO_SINGLE_LETTER_CODES = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
SINGLE_TO_THREE_LETTER_CODES = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}
ATOMIC_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "P": 30.974,
    "S": 32.065,
}
UNIPROT_ACCESSION_PATTERN = (
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}"
)
# %%
CHEMICAL_ELEMENT_PATTERN = r"""
    A[cglmrstu]     | 
    B[aehikr]?      | 
    C[adeflmnorsu]? | 
    D[bsy]          | 
    E[rsu]          | 
    F[elmr]?        | 
    G[ade]          | 
    H[efgos]?       | 
    I[nr]?          | 
    Kr?             | 
    L[airuv]        | 
    M[cdgnot]       | 
    N[abdehiop]?    | 
    O[gs]?          | 
    P[abdmortu]?    | 
    R[abefghnu]     | 
    S[bcegimnr]?    | 
    T[abcehilms]    | 
    U               | 
    V               | 
    W               | 
    Xe              | 
    Yb?             | 
    Z[nr]
""".replace("\n", "").replace(" ", "")
