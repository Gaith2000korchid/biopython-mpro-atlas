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


def _feature_span(feature: dict[str, Any]) -> tuple[int, int] | None:
    location = feature.get("location") or {}
    start = (location.get("start") or {}).get("value")
    end = (location.get("end") or {}).get("value")
    if start is None or end is None:
        return None
    start_i, end_i = int(start), int(end)
    if start_i > end_i:
        start_i, end_i = end_i, start_i
    return start_i, end_i


def domain_range_from_uniprot(entry: dict[str, Any], prosite_id: str = PROSITE_ID) -> dict[str, Any] | None:
    """Locate the M-pro interval on a UniProt entry from features, not from PDB lists."""
    for feature in entry.get("features") or []:
        for ref in feature.get("featureCrossReferences") or []:
            if ref.get("database") == "PROSITE" and ref.get("id") == prosite_id:
                span = _feature_span(feature)
                if span:
                    return {
                        "start": span[0],
                        "end": span[1],
                        "source": "uniprot_feature_prosite",
                        "feature_type": feature.get("type"),
                        "description": feature.get("description"),
                    }
    for feature in entry.get("features") or []:
        description = (feature.get("description") or "").lower()
        if any(keyword in description for keyword in _DOMAIN_KEYWORDS):
            span = _feature_span(feature)
            if span:
                return {
                    "start": span[0],
                    "end": span[1],
                    "source": "uniprot_feature_description",
                    "feature_type": feature.get("type"),
                    "description": feature.get("description"),
                }
    return None


def parse_chain_ranges(chains_value: str | None) -> list[dict[str, Any]]:
    """Parse UniProt PDB 'Chains' annotations such as 'A=3264-3569' or 'A/B=3264-3569'."""
    if not chains_value:
        return []
    assignments: list[dict[str, Any]] = []
    for match in _CHAIN_ASSIGNMENT.finditer(chains_value):
        start, end = int(match.group(2)), int(match.group(3))
        if start > end:
            start, end = end, start
        for chain_id in match.group(1).split("/"):
            chain_id = chain_id.strip()
            if chain_id:
                assignments.append({"chain": chain_id, "start": start, "end": end})
    return assignments


def inclusive_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> int:
    start = max(start_a, start_b)
    end = min(end_a, end_b)
    return max(0, end - start + 1)


def parse_resolution_angstrom(value: str | None) -> float | None:
    if not value:
        return None
    match = re.search(r"(\d+(?:\.\d+)?)", value)
    return float(match.group(1)) if match else None
