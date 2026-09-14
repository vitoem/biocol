from __future__ import annotations

import logging
from collections import Counter
from collections.abc import Iterable

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from biocol.sequence.alphabets import (
    GAP_LETTERS,
    NUCLEOTIDE_LETTERS,
    PROTEIN_ONLY_LETTERS,
    UNAMBIGUOUS_NUCLEOTIDE_LETTERS,
)

logger = logging.getLogger(__name__)

SequenceType = str
SequenceLike = str | Seq | SeqRecord


def _residues(sequence: SequenceLike) -> str:
    if isinstance(sequence, SeqRecord):
        text = str(sequence.seq)
    else:
        text = str(sequence)
    return "".join(char for char in text.upper() if char not in GAP_LETTERS)


def _sequence_id(sequence: SequenceLike) -> str:
    if isinstance(sequence, SeqRecord) and sequence.id:
        return sequence.id
    return "<no id>"


def detect_sequence_type(sequence: SequenceLike) -> SequenceType:
    """Classify a sequence as ``nucleotide`` or ``protein``.
    Uses ``Bio.Data.IUPACData``.
    """
    residues = _residues(sequence)
    if not residues:
        raise ValueError("Sequence has no classifiable residues")

    if any(char in PROTEIN_ONLY_LETTERS for char in residues):
        sequence_type = "protein"
    elif all(char in NUCLEOTIDE_LETTERS for char in residues) and any(
        char in UNAMBIGUOUS_NUCLEOTIDE_LETTERS for char in residues
    ):
        sequence_type = "nucleotide"
    else:
        sequence_type = "protein"

    logger.debug(
        "detect_sequence_type: id=%s length=%s type=%s",
        _sequence_id(sequence),
        len(residues),
        sequence_type,
    )
    return sequence_type


def detect_query_type(records: Iterable[SequenceLike]) -> SequenceType:
    """Type of a full FASTA: the majority type.

    Mixed classifications (e.g. short peptides that look like DNA) do not
    abort; protein wins ties.
    """
    records_list = list(records)
    types = [detect_sequence_type(record) for record in records_list]
    counts = Counter(types)
    protein_n = counts.get("protein", 0)
    nucleotide_n = counts.get("nucleotide", 0)
    query_type = "protein" if protein_n >= nucleotide_n else "nucleotide"
    if len(counts) > 1:
        minority_ids = [
            _sequence_id(record)
            for record, seq_type in zip(records_list, types, strict=True)
            if seq_type != query_type
        ]
        logger.debug(
            "detect_query_type: mixed FASTA; using %s (protein=%s nucleotide=%s) minority ids=%s",
            query_type,
            protein_n,
            nucleotide_n,
            minority_ids[:20],
        )
    logger.info(
        "Query/database FASTA classified as %s (%s sequences)",
        query_type,
        len(types),
    )
    return query_type
