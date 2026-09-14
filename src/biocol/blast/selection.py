from __future__ import annotations

import logging

from biocol.exceptions import BlastError

logger = logging.getLogger(__name__)

_DIRECT_PROGRAMS = {
    ("nucleotide", "protein"): "blastx",
    ("protein", "protein"): "blastp",
    ("protein", "nucleotide"): "tblastn",
}


def select_blast_program(
    query_type: str,
    database_type: str,
    *,
    translated: bool = False,
) -> str:
    """Choose blastn, blastp, blastx, tblastn, or tblastx.

    If query and database are nucleotide, the default is blastn.
    ``translated=True`` (explicit option, off by default) selects tblastx.
    In other combinations ``translated`` is ignored.
    """
    if query_type not in {"nucleotide", "protein"}:
        raise BlastError(f"Unsupported query type: {query_type}")
    if database_type not in {"nucleotide", "protein"}:
        raise BlastError(f"Unsupported database type: {database_type}")

    both_nucleotide = query_type == "nucleotide" and database_type == "nucleotide"
    if both_nucleotide:
        program = "tblastx" if translated else "blastn"
    else:
        program = _DIRECT_PROGRAMS[(query_type, database_type)]

    logger.debug(
        "select_blast_program: query=%s db=%s translated=%s → %s",
        query_type,
        database_type,
        translated,
        program,
    )
    return program
