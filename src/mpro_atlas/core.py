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
BASES = {
    "prosite": "https://prosite.expasy.org",
    "uniprot": "https://rest.uniprot.org/uniprotkb",
    "pdb": "https://files.rcsb.org/download",
}


def session_with_retries() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.8, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers.update({"User-Agent": "biopython-mpro-atlas/0.1 (educational research)"})
    return session


def download(session: requests.Session, url: str) -> str:
    response = session.get(url, timeout=30)
    response.raise_for_status()
    return response.text


def parse_aux(text: str) -> dict[str, Any]:
    """Extract profile-wide PDB IDs and matched UniProt accessions from PROSITE.AUX."""
    accession = re.search(r"^AC\s+(PS\d+);", text, re.MULTILINE)
    release = re.search(r"^NR\s+/RELEASE=([^;,]+)", text, re.MULTILINE)
    if accession is None:
        raise ValueError("Missing PROSITE accession in auxiliary record")
    pdb_ids: set[str] = set()
    matches: set[str] = set()
    for line in text.splitlines():
        if line.startswith("3D   "):
            pdb_ids.update(re.findall(r"\b[A-Za-z0-9]{4}\b(?=;)", line[5:]))
        if line.startswith("DR   "):
            matches.update(re.findall(r"\b([A-Z0-9]{6,10})\s*,\s*[^,;]+\s*,\s*[TNPF]\s*;", line[5:]))
    return {"accession": accession.group(1), "release": release.group(1) if release else None,
            "pdb_ids": sorted(pdb_ids), "uniprot_accessions": sorted(matches)}


def parse_uniprot(entry: dict[str, Any]) -> dict[str, Any]:
    accession = entry.get("primaryAccession")
    sequence = entry.get("sequence", {}).get("value")
    if not accession or not sequence:
        raise ValueError("UniProt entry lacks accession or sequence")
    cross_refs = entry.get("uniProtKBCrossReferences") or []
    pdb = {}
    prosite = set()
    for ref in cross_refs:
        if ref.get("database") == "PDB" and ref.get("id"):
            pdb[ref["id"].upper()] = {p["key"]: p["value"] for p in ref.get("properties", []) if "key" in p and "value" in p}
        if ref.get("database") == "PROSITE" and ref.get("id"):
            prosite.add(ref["id"])
    return {"accession": accession, "sequence": sequence,
            "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
            "organism": entry.get("organism", {}).get("scientificName"),
            "pdb": pdb, "prosite": sorted(prosite)}


def parse_cif(text: str, pdb_id: str) -> dict[str, Any]:
    structure = MMCIFParser(QUIET=True).get_structure(pdb_id, StringIO(text))
    model = next(structure.get_models())
    chains = list(model.get_chains())
    return {"pdb_id": pdb_id, "chain_count": len(chains),
            "chains": [{"id": chain.id, "residue_count": sum(1 for residue in chain if residue.id[0] == " ")}
                       for chain in chains]}


def build_report(session: requests.Session, accession: str = DEFAULT_ACCESSION, pdb_id: str | None = None,
                 include_structure: bool = True) -> tuple[dict[str, Any], SeqRecord]:
    if not re.fullmatch(r"[A-Z0-9]{6,10}", accession):
        raise ValueError("Invalid UniProt accession")
    primary = Prosite.read(StringIO(download(session, f"{BASES['prosite']}/{PROSITE_ID}.txt")))
    aux = parse_aux(download(session, f"{BASES['prosite']}/{PROSITE_ID}.aux"))
    if aux["accession"] != PROSITE_ID:
        raise ValueError("Auxiliary accession differs from requested profile")
    response = session.get(f"{BASES['uniprot']}/{accession}.json", timeout=30)
    response.raise_for_status()
    protein = parse_uniprot(response.json())
    if protein["accession"] != accession:
        raise ValueError("UniProt returned a different accession")
    shared = sorted(set(aux["pdb_ids"]) & protein["pdb"].keys())
    selected = pdb_id.upper() if pdb_id else (shared[0] if shared else None)
    if selected and selected not in protein["pdb"]:
        raise ValueError(f"{selected} is not cross-referenced by UniProt {accession}")
    report: dict[str, Any] = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "profile": {"id": PROSITE_ID, "description": primary.description,
                    "primary_record_pdb_count": len(primary.pdb_structs),
                    "aux_release": aux["release"], "aux_pdb_count": len(aux["pdb_ids"]),
                    "uniprot_accession_in_aux": accession in aux["uniprot_accessions"]},
        "protein": {k: v for k, v in protein.items() if k != "sequence"},
        "shared_pdb_ids": shared,
        "selected_pdb": selected,
        "selected_pdb_metadata": protein["pdb"].get(selected) if selected else None,
        "caveat": "PROSITE.AUX 3D entries are profile-wide; overlap with UniProt PDB cross-references is a database cross-reference, not a residue-level mapping or an independent validation of a ScanProsite match.",
        "sources": {"prosite_primary": f"{BASES['prosite']}/{PROSITE_ID}.txt",
                    "prosite_aux": f"{BASES['prosite']}/{PROSITE_ID}.aux",
                    "uniprot": f"{BASES['uniprot']}/{accession}.json"},
    }
    if selected and include_structure:
        url = f"{BASES['pdb']}/{selected}.cif"
        report["structure"] = parse_cif(download(session, url), selected)
        report["sources"]["cif"] = url
    fasta = SeqRecord(Seq(protein["sequence"]), id=accession,
                      description=f"{protein['protein_name'] or 'protein'} | {protein['organism'] or 'unknown organism'}")
    return report, fasta
