"""Custom MSATK exceptions with user-facing messages."""

from __future__ import annotations


class MSATKError(Exception):
    """Base class for MSATK errors."""


class AlignmentFormatError(MSATKError):
    """Raised when MSATK cannot parse or detect an alignment format."""


class AlignmentValidationError(MSATKError):
    """Raised when an alignment fails strict validation."""


class CodonAlignmentError(MSATKError):
    """Raised when codon-aware analysis cannot be performed safely."""
