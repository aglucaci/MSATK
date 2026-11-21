"""Alignment data models, unified loading and the format capability registry.

>>> from msatk import load_alignment
>>> aln = load_alignment("genomes.maf.gz")
>>> aln.kind, aln.summary()["blocks"]
"""

from .blocks import AlignmentBlock, BlockAlignment, BlockRow
from .graphs import GraphAlignment
from .loader import SequenceAlignment, detect_input_format, inspect_alignment, load_alignment
from .mappings import MappingAlignment, MappingRecord
from .profiles import ProfileModel
from .reads import ReadAlignment, ReadRecord
from .registry import FORMATS, KINDS, capability_registry, compatibility_markdown, get_format

__all__ = [
    "AlignmentBlock",
    "BlockAlignment",
    "BlockRow",
    "FORMATS",
    "GraphAlignment",
    "KINDS",
    "MappingAlignment",
    "MappingRecord",
    "ProfileModel",
    "ReadAlignment",
    "ReadRecord",
    "SequenceAlignment",
    "capability_registry",
    "compatibility_markdown",
    "detect_input_format",
    "get_format",
    "inspect_alignment",
    "load_alignment",
]
