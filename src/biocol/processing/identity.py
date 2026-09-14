"""Full-query identity from BLAST HSPs using aligned ``qseq`` / ``sseq``."""

from __future__ import annotations

from collections.abc import Iterable, Mapping

from biocol.sequence.alphabets import GAP_LETTERS

_GAP = GAP_LETTERS | {"*"}


def _query_step(qseq: str, qstart: int, qend: int, query_is_nucleotide: bool) -> int:
    """Nucleotides per aligned query residue (3 for blastx/tblastx)."""
    if not query_is_nucleotide:
        return 1
    residues = sum(1 for char in qseq if char.upper() not in GAP_LETTERS)
    if residues == 0:
        return 1
    span = abs(qend - qstart) + 1
    if span == 3 * residues:
        return 3
    return 1


def identical_query_positions(
    qseq: str,
    sseq: str,
    qstart: int,
    qend: int,
    *,
    query_is_nucleotide: bool,
) -> set[int]:
    """1-based query coordinates that match the subject in this HSP."""
    if len(qseq) != len(sseq):
        return set()
    step = _query_step(qseq, qstart, qend, query_is_nucleotide)
    plus = qstart <= qend
    position = qstart
    matched: set[int] = set()
    for query_char, subject_char in zip(qseq, sseq, strict=True):
        query_u = query_char.upper()
        subject_u = subject_char.upper()
        if query_u in GAP_LETTERS:
            continue
        if query_u == subject_u and query_u not in _GAP and subject_u not in GAP_LETTERS:
            if plus:
                matched.update(range(position, position + step))
            else:
                matched.update(range(position - step + 1, position + 1))
        if plus:
            position += step
        else:
            position -= step
    return matched


def full_query_identity_percent(
    hsps: Iterable[Mapping[str, object]],
    query_length: int,
    *,
    query_is_nucleotide: bool,
) -> float | None:
    """Union of identical query positions across HSPs, divided by query length.

    ``hsps`` must belong to the same query–subject pair. Returns ``None``
    when the value cannot be computed (no length, no aligned sequences).
    """
    if query_length <= 0:
        return None
    union: set[int] = set()
    saw_alignment = False
    for hsp in hsps:
        qseq = hsp.get("qseq")
        sseq = hsp.get("sseq")
        qstart = hsp.get("qstart")
        qend = hsp.get("qend")
        if qseq is None or sseq is None or qstart is None or qend is None:
            continue
        if str(qseq) in {"", "nan", "<NA>"} or str(sseq) in {"", "nan", "<NA>"}:
            continue
        try:
            start = int(float(str(qstart)))
            end = int(float(str(qend)))
        except (TypeError, ValueError):
            continue
        saw_alignment = True
        union |= identical_query_positions(
            str(qseq),
            str(sseq),
            start,
            end,
            query_is_nucleotide=query_is_nucleotide,
        )
    if not saw_alignment:
        return None
    clipped = {pos for pos in union if 1 <= pos <= query_length}
    return round(100.0 * len(clipped) / query_length, 2)
