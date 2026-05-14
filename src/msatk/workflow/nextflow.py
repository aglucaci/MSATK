"""Nextflow integration helpers."""


def nextflow_process() -> str:
    return """process MSATK_PROFILE {
  input:
  path alignment

  output:
  path "msatk_results"

  script:
  '''
  msatk profile ${alignment} --out msatk_results --force
  '''
}
"""
