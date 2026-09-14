from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from biocol.exceptions import EmptyFastaError, InvalidFastaError
from biocol.sequence.alphabets import VALID_RESIDUES

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna", ".faa", ".fas"}


def check_fasta_path(path: str | Path) -> Path:
    """Check that the path exists, is a file, and has a FASTA extension."""
    fasta_path = Path(path)
    if not fasta_path.exists():
        raise FileNotFoundError(f"FASTA file not found: {fasta_path}")
    if not fasta_path.is_file():
        raise InvalidFastaError(f"Path is not a file: {fasta_path}")
    if fasta_path.suffix.lower() not in FASTA_EXTENSIONS:
        raise InvalidFastaError(
            f"Unsupported extension: {fasta_path.suffix}. "
            f"Use one of: {', '.join(sorted(FASTA_EXTENSIONS))}"
        )
    return fasta_path


def validate_seq_record(record: SeqRecord, source: str | Path) -> None:
    """Check that a Biopython SeqRecord has a usable sequence."""
    sequence = str(record.seq).upper()
    if not sequence:
        raise InvalidFastaError(f"Sequence '{record.id}' is empty in {source}")
    invalid = {char for char in sequence if char not in VALID_RESIDUES}
    if invalid:
        raise InvalidFastaError(
            f"Sequence '{record.id}' contains invalid characters: "
            f"{''.join(sorted(invalid))}"
        )


def validate_fasta_file(path: str | Path) -> Path:
    """Check existence, extension, and valid FASTA content.

    Parses the file once with ``Bio.SeqIO``.
    """
    from biocol.sequence.reader import read_fasta

    fasta_path = check_fasta_path(path)
    read_fasta(fasta_path)
    return fasta_path
