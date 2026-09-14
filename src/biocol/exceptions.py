class FastaError(ValueError):
    """Base error for FASTA files."""


class EmptyFastaError(FastaError):
    """The file exists but contains no sequences."""


class InvalidFastaError(FastaError):
    """The file is not a valid FASTA or the sequences are not acceptable."""


class MixedSequenceTypeError(FastaError):
    """A multifasta mixes nucleotide and protein sequences."""


class BlastError(ValueError):
    """Base error for BLAST selection or execution."""


class DatabaseError(BlastError):
    """The FASTA database does not exist, is unrecognized, or is invalid."""


class MixedDatabaseTypeError(DatabaseError):
    """A folder or list mixes nucleotide and protein databases."""


class BlastExecutionError(BlastError):
    """makeblastdb or BLAST+ execution failed."""


class MetadataError(ValueError):
    """The accessions/descriptors file is invalid."""


class HmmError(ValueError):
    """Base error for HMMER / hmmscan."""


class HmmExecutionError(HmmError):
    """hmmpress or hmmscan execution failed."""


class DiamondError(ValueError):
    """Base error for Diamond search."""


class DiamondExecutionError(DiamondError):
    """diamond makedb or blastp execution failed."""


class KofamError(ValueError):
    """Base error for KofamScan."""


class KofamExecutionError(KofamError):
    """exec_annotation execution failed."""
