"""Workflow integration snippets."""

from .nextflow import nextflow_process
from .snakemake import snakemake_rule

__all__ = ["snakemake_rule", "nextflow_process"]
