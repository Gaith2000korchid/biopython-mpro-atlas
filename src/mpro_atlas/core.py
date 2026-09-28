"""Network-bound retrieval separated from deterministic parsing and analysis."""

from __future__ import annotations

from datetime import datetime, timezone
from io import StringIO
import re
from typing import Any

from Bio.ExPASy import Prosite
from Bio.PDB import MMCIFParser
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry


PROSITE_ID = "PS51442"  # Coronavirus main protease domain profile
DEFAULT_ACCESSION = "P0DTD1"  # SARS-CoV-2 replicase polyprotein 1ab
# Inclusive polyprotein coordinates of mature M-pro / nsp5 on P0DTD1.
# Used only when UniProt features do not locate PS51442.
DEFAULT_DOMAIN_RANGES = {
    DEFAULT_ACCESSION: (3264, 3569),
}
BASES = {
    "prosite": "https://prosite.expasy.org",
    "uniprot": "https://rest.uniprot.org/uniprotkb",
    "pdb": "https://files.rcsb.org/download",
}

_PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_ACCESSION = re.compile(r"^[A-Z0-9]{6,10}$")
_AUX_AC = re.compile(r"^AC\s+(PS\d+);", re.MULTILINE)
_AUX_RELEASE = re.compile(r"^NR\s+/RELEASE=([^;,]+)", re.MULTILINE)
_AUX_DR = re.compile(r"\b([A-Z0-9]{6,10})\s*,\s*[^,;]+\s*,\s*[TNPF]\s*;")
_CHAIN_ASSIGNMENT = re.compile(
    r"([A-Za-z0-9][A-Za-z0-9/]*)\s*=\s*(\d+)\s*-\s*(\d+)"
)
_DOMAIN_KEYWORDS = (
    "3c-like",
    "3cl-like",
    "3clpro",
    "main protease",
    "m-pro",
    "mpro",
    "nsp5",
)


def session_with_retries() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.8, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers.update({"User-Agent": "biopython-mpro-atlas/0.2 (educational research)"})
    return session


def download(session: requests.Session, url: str) -> str:
    response = session.get(url, timeout=30)
    response.raise_for_status()
    return response.text


def parse_aux(text: str) -> dict[str, Any]:
    """Extract profile-wide PDB IDs and matched UniProt accessions from PROSITE.AUX."""
    if not text or not text.strip():
        raise ValueError("Empty PROSITE auxiliary record")
    accession = _AUX_AC.search(text)
    release = _AUX_RELEASE.search(text)
    if accession is None:
        raise ValueError("Missing PROSITE accession in auxiliary record")

    pdb_ids: set[str] = set()
    matches: set[str] = set()
    saw_3d = False
    saw_dr = False
    for line in text.splitlines():
        if line.startswith("3D   "):
            saw_3d = True
            for token in line[5:].split(";"):
                token = token.strip()
                if _PDB_ID.fullmatch(token):
                    pdb_ids.add(token.upper())
        elif line.startswith("DR   "):
            saw_dr = True
            matches.update(_AUX_DR.findall(line[5:]))

    if not saw_3d and not saw_dr:
        raise ValueError("Auxiliary record has no 3D or DR lines")

    return {
        "accession": accession.group(1),
        "release": release.group(1) if release else None,
        "pdb_ids": sorted(pdb_ids),
        "uniprot_accessions": sorted(matches),
    }
