"""Network-bound retrieval separated from deterministic parsing and analysis."""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from importlib.metadata import PackageNotFoundError, version as pkg_version
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

from .mapping import (
    attach_sifts,
    coverage_counts,
    fetch_sifts,
    filter_domain_structures,
    inclusive_overlap,
    map_pdb_to_domain,
    parse_chain_ranges,
    parse_resolution_angstrom,
    parse_sifts_uniprot,
    select_structure,
    sifts_domain_coverage,
)

PROSITE_ID = "PS51442"
DEFAULT_ACCESSION = "P0DTD1"
DEFAULT_DOMAIN_RANGES = {DEFAULT_ACCESSION: (3264, 3569)}
DEFAULT_MIN_COVERAGE = 50.0
PACKAGE_VERSION = "0.3.0"
BASES = {
    "prosite": "https://prosite.expasy.org",
    "uniprot": "https://rest.uniprot.org/uniprotkb",
    "pdb": "https://files.rcsb.org/download",
    "pdbe": "https://www.ebi.ac.uk/pdbe/api",
}
_PDB_ID = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
_ACCESSION = re.compile(r"^[A-Z0-9]{6,10}$")
_AUX_AC = re.compile(r"^AC\s+(PS\d+);", re.MULTILINE)
_AUX_RELEASE = re.compile(r"^NR\s+/RELEASE=([^;,]+)", re.MULTILINE)
_AUX_DR = re.compile(r"\b([A-Z0-9]{6,10})\s*,\s*[^,;]+\s*,\s*[TNPF]\s*;")
_DOMAIN_KEYWORDS = ("3c-like", "3cl-like", "3clpro", "main protease", "m-pro", "mpro", "nsp5")


def session_with_retries() -> requests.Session:
    session = requests.Session()
    retries = Retry(total=3, backoff_factor=0.8, status_forcelist=[429, 500, 502, 503, 504])
    session.mount("https://", HTTPAdapter(max_retries=retries))
    session.headers.update({"User-Agent": f"biopython-mpro-atlas/{PACKAGE_VERSION} (educational research)"})
    return session


def download(session: requests.Session, url: str) -> str:
    response = session.get(url, timeout=30)
    response.raise_for_status()
    return response.text


def parse_aux(text: str) -> dict[str, Any]:
    if not text or not text.strip():
        raise ValueError("Empty PROSITE auxiliary record")
    accession = _AUX_AC.search(text)
    release = _AUX_RELEASE.search(text)
    if accession is None:
        raise ValueError("Missing PROSITE accession in auxiliary record")
    pdb_ids: set[str] = set()
    matches: set[str] = set()
    saw_3d = saw_dr = False
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
    return {"accession": accession.group(1), "release": release.group(1) if release else None,
            "pdb_ids": sorted(pdb_ids), "uniprot_accessions": sorted(matches)}


def _feature_span(feature: dict[str, Any]) -> tuple[int, int] | None:
    location = feature.get("location") or {}
    start = (location.get("start") or {}).get("value")
    end = (location.get("end") or {}).get("value")
    if start is None or end is None:
        return None
    start_i, end_i = int(start), int(end)
    return (end_i, start_i) if start_i > end_i else (start_i, end_i)


def domain_range_from_uniprot(entry: dict[str, Any], prosite_id: str = PROSITE_ID) -> dict[str, Any] | None:
    for feature in entry.get("features") or []:
        for ref in feature.get("featureCrossReferences") or []:
            if ref.get("database") == "PROSITE" and ref.get("id") == prosite_id:
                span = _feature_span(feature)
                if span:
                    return {"start": span[0], "end": span[1], "source": "uniprot_feature_prosite",
                            "feature_type": feature.get("type"), "description": feature.get("description")}
    for feature in entry.get("features") or []:
        description = (feature.get("description") or "").lower()
        if any(keyword in description for keyword in _DOMAIN_KEYWORDS):
            span = _feature_span(feature)
            if span:
                return {"start": span[0], "end": span[1], "source": "uniprot_feature_description",
                        "feature_type": feature.get("type"), "description": feature.get("description")}
    return None


def parse_uniprot(entry: dict[str, Any]) -> dict[str, Any]:
    accession = entry.get("primaryAccession")
    sequence = entry.get("sequence", {}).get("value")
    if not accession or not sequence:
        raise ValueError("UniProt entry lacks accession or sequence")
    pdb, prosite = {}, set()
    for ref in entry.get("uniProtKBCrossReferences") or []:
        if ref.get("database") == "PDB" and ref.get("id"):
            pdb[ref["id"].upper()] = {p["key"]: p["value"] for p in ref.get("properties", []) if "key" in p and "value" in p}
        if ref.get("database") == "PROSITE" and ref.get("id"):
            prosite.add(ref["id"])
    return {"accession": accession, "sequence": sequence,
            "protein_name": entry.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value"),
            "organism": entry.get("organism", {}).get("scientificName"),
            "pdb": pdb, "prosite": sorted(prosite), "domain": domain_range_from_uniprot(entry)}


def parse_cif(text: str, pdb_id: str) -> dict[str, Any]:
    structure = MMCIFParser(QUIET=True).get_structure(pdb_id, StringIO(text))
    model = next(structure.get_models())
    chains = list(model.get_chains())
    return {"pdb_id": pdb_id, "chain_count": len(chains),
            "chains": [{"id": chain.id, "residue_count": sum(1 for residue in chain if residue.id[0] == " ")} for chain in chains]}


def resolve_domain(protein, domain_start=None, domain_end=None):
    if domain_start is not None and domain_end is not None:
        start, end = (domain_end, domain_start) if domain_start > domain_end else (domain_start, domain_end)
        return {"start": start, "end": end, "source": "cli", "feature_type": None, "description": "User-supplied domain coordinates"}
    if protein.get("domain"):
        return protein["domain"]
    fallback = DEFAULT_DOMAIN_RANGES.get(protein["accession"])
    if fallback:
        return {"start": fallback[0], "end": fallback[1], "source": "accession_default", "feature_type": None,
                "description": "Documented nsp5 / M-pro interval on P0DTD1"}
    return None


def _csv_cell(value: Any) -> str:
    text = "" if value is None else str(value)
    if any(char in text for char in [",", '"', "\n"]):
        return '"' + text.replace('"', '""') + '"'
    return text


def domain_structures_csv(rows):
    header = ("pdb_id,method,resolution,chain,uniprot_start,uniprot_end,uniprot_overlap,"
              "uniprot_coverage_pct,sifts_overlap,sifts_coverage_pct,sifts_missing,ligands")
    lines = [header]
    for row in rows:
        sifts = row.get("sifts") or {}
        targets = row["overlapping_chains"] or [{"chain": "", "start": "", "end": "", "overlap_residues": row["overlap_residues"]}]
        for chain in targets:
            lines.append(",".join([
                row["pdb_id"], _csv_cell(row.get("method")), _csv_cell(row.get("resolution")),
                _csv_cell(chain.get("chain")), str(chain.get("start", "")), str(chain.get("end", "")),
                str(chain.get("overlap_residues", "")), str(row.get("coverage_pct", "")),
                str(sifts.get("overlap_residues", "")), str(sifts.get("coverage_pct", "")),
                str(sifts.get("missing_residues", "")), _csv_cell(row.get("ligands")),
            ]))
    return "\n".join(lines) + "\n"


def _digest(text: str) -> str:
    return sha256(text.encode("utf-8")).hexdigest()


def _biopython_version() -> str | None:
    try:
        return pkg_version("biopython")
    except PackageNotFoundError:
        return getattr(__import__("Bio"), "__version__", None)


def build_report(session, accession=DEFAULT_ACCESSION, pdb_id=None, include_structure=True,
                 domain_start=None, domain_end=None, min_coverage=DEFAULT_MIN_COVERAGE,
                 include_sifts=True, sifts_limit=1):
    if not _ACCESSION.fullmatch(accession):
        raise ValueError("Invalid UniProt accession")
    if min_coverage < 0 or min_coverage > 100:
        raise ValueError("min_coverage must be between 0 and 100")
    if (domain_start is None) ^ (domain_end is None):
        raise ValueError("domain_start and domain_end must be supplied together")
    if sifts_limit < 0:
        raise ValueError("sifts_limit must be >= 0")

    primary_text = download(session, f"{BASES['prosite']}/{PROSITE_ID}.txt")
    aux_text = download(session, f"{BASES['prosite']}/{PROSITE_ID}.aux")
    primary = Prosite.read(StringIO(primary_text))
    aux = parse_aux(aux_text)
    if aux["accession"] != PROSITE_ID:
        raise ValueError("Auxiliary accession differs from requested profile")
    uniprot_url = f"{BASES['uniprot']}/{accession}.json"
    response = session.get(uniprot_url, timeout=30)
    response.raise_for_status()
    uniprot_text = response.text
    protein = parse_uniprot(response.json())
    if protein["accession"] != accession:
        raise ValueError("UniProt returned a different accession")

    shared = sorted(set(aux["pdb_ids"]) & protein["pdb"].keys())
    domain = resolve_domain(protein, domain_start, domain_end)
    mapped, unmapped = [], 0
    if domain:
        mapped = [map_pdb_to_domain(pdb, protein["pdb"][pdb], domain["start"], domain["end"]) for pdb in shared]
        unmapped = sum(1 for row in mapped if row["unmapped"])
    tiers = coverage_counts(mapped)
    domain_rows = [row for row in mapped if row["covers_domain"] and row["coverage_pct"] >= min_coverage]
    domain_rows.sort(key=lambda row: (-row["coverage_pct"], row["resolution_angstrom"] if row["resolution_angstrom"] is not None else 1e9, row["pdb_id"]))
    selected, selected_reason = select_structure(pdb_id, protein, domain, domain_rows, min_coverage)

    sifts_urls = []
    if include_sifts and domain and sifts_limit:
        targets = []
        if selected:
            targets.append(selected)
        for row in domain_rows:
            if row["pdb_id"] not in targets:
                targets.append(row["pdb_id"])
            if len(targets) >= sifts_limit:
                break
        enriched = {row["pdb_id"]: row for row in domain_rows}
        for pdb in targets:
            payload, url = fetch_sifts(session, pdb, BASES["pdbe"])
            sifts_urls.append(url)
            base = enriched.get(pdb) or map_pdb_to_domain(pdb, protein["pdb"][pdb], domain["start"], domain["end"])
            attached = attach_sifts(base, payload, accession, domain["start"], domain["end"])
            if pdb in enriched:
                enriched[pdb] = attached
        domain_rows = [enriched[row["pdb_id"]] for row in domain_rows]

    report = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "profile": {"id": PROSITE_ID, "description": primary.description, "primary_record_pdb_count": len(primary.pdb_structs),
                    "aux_release": aux["release"], "aux_pdb_count": len(aux["pdb_ids"]),
                    "uniprot_accession_in_aux": accession in aux["uniprot_accessions"]},
        "protein": {k: v for k, v in protein.items() if k != "sequence"},
        "domain": domain, "min_coverage": min_coverage, "coverage_tiers": tiers,
        "shared_pdb_ids": shared, "shared_pdb_count": len(shared),
        "domain_structure_count": len(domain_rows), "unmapped_shared_pdb_count": unmapped,
        "domain_structures": domain_rows, "selected_pdb": selected,
        "selected_pdb_metadata": protein["pdb"].get(selected) if selected else None,
        "selected_reason": selected_reason,
        "caveat": "UniProt Chains comparison is an annotation-range test. SIFTS/PDBe segments map PDB polymer ranges to UniProt positions and can report missing domain residues, but they are still segment alignments, not a per-atom experimental proof of activity.",
        "provenance": {
            "mpro_atlas_version": PACKAGE_VERSION,
            "biopython_version": _biopython_version(),
            "prosite_release": aux["release"],
            "prosite_primary_sha256": _digest(primary_text),
            "prosite_aux_sha256": _digest(aux_text),
            "uniprot_sha256": _digest(uniprot_text),
        },
        "sources": {"prosite_primary": f"{BASES['prosite']}/{PROSITE_ID}.txt",
                    "prosite_aux": f"{BASES['prosite']}/{PROSITE_ID}.aux", "uniprot": uniprot_url},
    }
    if sifts_urls:
        report["sources"]["pdbe_sifts"] = sifts_urls
    if selected and include_structure:
        url = f"{BASES['pdb']}/{selected}.cif"
        report["structure"] = parse_cif(download(session, url), selected)
        report["sources"]["cif"] = url
    fasta = SeqRecord(Seq(protein["sequence"]), id=accession,
                      description=f"{protein['protein_name'] or 'protein'} | {protein['organism'] or 'unknown organism'}")
    return report, fasta, domain_structures_csv(domain_rows)
