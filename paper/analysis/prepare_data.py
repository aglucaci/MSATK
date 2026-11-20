"""Download and prepare the public alignments used in the MSATK Application Note.

Datasets
--------
* Pfam seed alignments (InterPro API): PF00042 globin, PF00069 protein kinase,
  PF00071 Ras family, PF00076 RRM. Stockholm is converted to FASTA; Pfam insert
  columns ('.') are written as gaps ('-') and residues upper-cased.
* Influenza A(H3N2) haemagglutinin coding sequences (NCBI nuccore, first 400
  records of the query below, CDS extracted with ``rettype=fasta_cds_na``). All
  retrieved CDS are 1,701 nt, so they form an ungapped codon alignment without
  re-alignment. Duplicate sequences are deliberately retained so that QC can
  report them.

Usage::

    python paper/analysis/prepare_data.py            # download + convert
    python paper/analysis/prepare_data.py --offline  # convert existing files only
"""

from __future__ import annotations

import argparse
import gzip
import re
import urllib.parse
import urllib.request
from pathlib import Path

from Bio import AlignIO, SeqIO

DATA = Path(__file__).parent / "data"
PFAM = ["PF00042", "PF00069", "PF00071", "PF00076"]
NCBI_TERM = "Influenza A virus[orgn] AND H3N2 AND hemagglutinin AND complete cds AND 2017[pdat]"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def _get(url: str, data: dict[str, str] | None = None) -> bytes:
    body = urllib.parse.urlencode(data).encode() if data else None
    with urllib.request.urlopen(url, data=body, timeout=120) as response:
        return response.read()


def download() -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    for acc in PFAM:
        url = f"https://www.ebi.ac.uk/interpro/wwwapi//entry/pfam/{acc}/?annotation=alignment:seed"
        (DATA / f"{acc}.seed.sto.gz").write_bytes(_get(url))
    query = urllib.parse.urlencode(
        {"db": "nuccore", "term": NCBI_TERM, "retmax": "400", "tool": "msatk_paper"}
    )
    ids = re.findall(r"<Id>(\d+)</Id>", _get(f"{EUTILS}/esearch.fcgi?{query}").decode())
    cds = _get(
        f"{EUTILS}/efetch.fcgi",
        {
            "db": "nuccore",
            "rettype": "fasta_cds_na",
            "retmode": "text",
            "id": ",".join(ids),
            "tool": "msatk_paper",
        },
    )
    (DATA / "h3n2_cds_raw.fasta").write_bytes(cds)


def convert() -> None:
    for acc in PFAM:
        with gzip.open(DATA / f"{acc}.seed.sto.gz", "rt") as handle:
            aln = AlignIO.read(handle, "stockholm")
        with (DATA / f"{acc}.fasta").open("w") as out:
            for rec in aln:
                seq = str(rec.seq).replace(".", "-").upper()
                out.write(f">{rec.id.replace('/', '_')}\n{seq}\n")
    with (DATA / "h3n2_ha.fasta").open("w") as out:
        for rec in SeqIO.parse(DATA / "h3n2_cds_raw.fasta", "fasta"):
            acc = rec.id.split("|")[1].split("_cds")[0]
            out.write(f">{acc}\n{str(rec.seq).upper()}\n")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("--offline", action="store_true")
    if not ap.parse_args().offline:
        download()
    convert()
