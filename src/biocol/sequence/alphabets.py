"""IUPAC alphabets from Biopython (``Bio.Data.IUPACData``)."""

from Bio.Data.IUPACData import (
    ambiguous_dna_letters,
    ambiguous_rna_letters,
    extended_protein_letters,
    unambiguous_dna_letters,
    unambiguous_rna_letters,
)

# DNA/RNA with ambiguity codes (K, R, Y, ...). These overlap amino acids.
NUCLEOTIDE_LETTERS = set((ambiguous_dna_letters + ambiguous_rna_letters).upper())
# Unambiguous bases: not confused with lysine (K), arginine (R), etc.
UNAMBIGUOUS_NUCLEOTIDE_LETTERS = set(
    (unambiguous_dna_letters + unambiguous_rna_letters).upper()
)
PROTEIN_LETTERS = set(extended_protein_letters.upper())
GAP_LETTERS = {"-", "."}
VALID_RESIDUES = NUCLEOTIDE_LETTERS | PROTEIN_LETTERS | GAP_LETTERS | {"*"}
# E, F, I, L, P, Q, X, O, Z, J, *
PROTEIN_ONLY_LETTERS = (PROTEIN_LETTERS - NUCLEOTIDE_LETTERS) | {"*"}
