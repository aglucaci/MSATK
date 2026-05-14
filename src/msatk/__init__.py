"""MSATK: Multiple Sequence Alignment Toolkit."""

from .api import profile_alignment
from .codon.profiler import CodonProfiler
from .core.profiler import MSATK, AlignmentProfiler
from .protein.profiler import ProteinProfiler
from .version import __version__

msatk = MSATK

__all__ = [
    "AlignmentProfiler",
    "MSATK",
    "msatk",
    "CodonProfiler",
    "ProteinProfiler",
    "profile_alignment",
    "__version__",
]
