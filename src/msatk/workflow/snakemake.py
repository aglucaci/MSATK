"""Snakemake integration helpers."""


def snakemake_rule() -> str:
    return """rule msatk_profile:
    input:
        "alignment.fasta"
    output:
        directory("alignment_msatk")
    shell:
        "msatk profile {input} --out {output} --force"
"""
