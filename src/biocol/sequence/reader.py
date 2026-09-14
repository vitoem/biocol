from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from biocol.exceptions import EmptyFastaError, InvalidFastaError
from biocol.sequence.validator import check_fasta_path, validate_seq_record


def read_fasta(path: str | Path) -> list[SeqRecord]:
    """Read a FASTA or multifasta with ``Bio.SeqIO.parse``.

    Returns Biopython ``SeqRecord`` objects (id, description, seq).
    """
    fasta_path = check_fasta_path(path)
    try:
        records = list(SeqIO.parse(fasta_path, "fasta"))
    except Exception as exc:
        raise InvalidFastaError(f"Could not parse FASTA: {fasta_path}") from exc

    if not records:
        raise EmptyFastaError(f"FASTA contains no sequences: {fasta_path}")

    for record in records:
        validate_seq_record(record, fasta_path)
        record.seq = record.seq.upper()
    return records
