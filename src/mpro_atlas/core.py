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
DEFAULT_DOMAIN_RANGES = {DEFAULT_ACCESSION: (3264, 3569)}
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
_CHAIN_ASSIGNMENT = re.compile(r"([A-Za-z0-9][A-Za-z0-9/]*)\s*=\s*(\d+)\s*-\s*(\d+)")
_DOMAIN_KEYWORDS = ("3c-like", "3cl-like", "3clpro", "main protease", "m-pro", "mpro", "nsp5")


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


def parse_chain_ranges(chains_value: str | None) -> list[dict[str, Any]]:
    if not chains_value:
        return []
    assignments = []
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
    return max(0, min(end_a, end_b) - max(start_a, start_b) + 1)


def parse_resolution_angstrom(value: str | None) -> float | None:
    if not value:
        return None
    match = re.search(r"(\d+(?:\.\d+)?)", value)
    return float(match.group(1)) if match else None


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


def map_pdb_to_domain(pdb_id: str, properties: dict[str, str], domain_start: int, domain_end: int) -> dict[str, Any]:
    domain_length = domain_end - domain_start + 1
    assignments = parse_chain_ranges(properties.get("Chains"))
    overlapping, best_overlap = [], 0
    for item in assignments:
        overlap = inclusive_overlap(item["start"], item["end"], domain_start, domain_end)
        if overlap:
            overlapping.append({**item, "overlap_residues": overlap})
            best_overlap = max(best_overlap, overlap)
    coverage = round(100.0 * best_overlap / domain_length, 2) if domain_length else 0.0
    return {"pdb_id": pdb_id, "method": properties.get("Method"), "resolution": properties.get("Resolution"),
            "resolution_angstrom": parse_resolution_angstrom(properties.get("Resolution")),
            "ligands": properties.get("Ligands"), "chains_annotation": properties.get("Chains"),
            "chain_ranges": assignments, "overlapping_chains": overlapping,
            "overlap_residues": best_overlap, "coverage_pct": coverage,
            "covers_domain": best_overlap > 0, "mapping_source": "uniprot_pdb_chains", "unmapped": not assignments}


def filter_domain_structures(pdb_annotations, domain_start, domain_end, candidate_ids=None, min_coverage=0.0):
    ids = candidate_ids if candidate_ids is not None else sorted(pdb_annotations)
    rows = []
    for pdb_id in ids:
        properties = pdb_annotations.get(pdb_id)
        if properties is None:
            continue
        row = map_pdb_to_domain(pdb_id, properties, domain_start, domain_end)
        if row["covers_domain"] and row["coverage_pct"] >= min_coverage:
            rows.append(row)
    rows.sort(key=lambda row: (-row["coverage_pct"], row["resolution_angstrom"] if row["resolution_angstrom"] is not None else 1e9, row["pdb_id"]))
    return rows


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
    header = "pdb_id,method,resolution,chain,uniprot_start,uniprot_end,overlap_residues,coverage_pct,ligands"
    lines = [header]
    for row in rows:
        targets = row["overlapping_chains"] or [{"chain": "", "start": "", "end": "", "overlap_residues": row["overlap_residues"]}]
        for chain in targets:
            lines.append(",".join([row["pdb_id"], _csv_cell(row.get("method")), _csv_cell(row.get("resolution")),
                                   _csv_cell(chain.get("chain")), str(chain.get("start", "")), str(chain.get("end", "")),
                                   str(chain.get("overlap_residues", "")), str(row.get("coverage_pct", "")),
                                   _csv_cell(row.get("ligands"))]))
    return "\n".join(lines) + "\n"


def build_report(session, accession=DEFAULT_ACCESSION, pdb_id=None, include_structure=True,
                 domain_start=None, domain_end=None, min_coverage=0.0):
    if not _ACCESSION.fullmatch(accession):
        raise ValueError("Invalid UniProt accession")
    if min_coverage < 0 or min_coverage > 100:
        raise ValueError("min_coverage must be between 0 and 100")
    if (domain_start is None) ^ (domain_end is None):
        raise ValueError("domain_start and domain_end must be supplied together")
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
    domain = resolve_domain(protein, domain_start, domain_end)
    domain_rows, unmapped = [], 0
    if domain:
        mapped = [map_pdb_to_domain(pdb, protein["pdb"][pdb], domain["start"], domain["end"]) for pdb in shared if pdb in protein["pdb"]]
        unmapped = sum(1 for row in mapped if row["unmapped"])
        domain_rows = [row for row in mapped if row["covers_domain"] and row["coverage_pct"] >= min_coverage]
        domain_rows.sort(key=lambda row: (-row["coverage_pct"], row["resolution_angstrom"] if row["resolution_angstrom"] is not None else 1e9, row["pdb_id"]))
    selected = None
    if pdb_id:
        selected = pdb_id.upper()
        if selected not in protein["pdb"]:
            raise ValueError(f"{selected} is not cross-referenced by UniProt {accession}")
        if domain_rows and selected not in {row["pdb_id"] for row in domain_rows}:
            raise ValueError(f"{selected} does not overlap the M-pro interval {domain['start']}-{domain['end']} at min_coverage={min_coverage}")
    elif domain_rows:
        selected = domain_rows[0]["pdb_id"]
    elif shared:
        selected = shared[0]
    report = {
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "profile": {"id": PROSITE_ID, "description": primary.description, "primary_record_pdb_count": len(primary.pdb_structs),
                    "aux_release": aux["release"], "aux_pdb_count": len(aux["pdb_ids"]),
                    "uniprot_accession_in_aux": accession in aux["uniprot_accessions"]},
        "protein": {k: v for k, v in protein.items() if k != "sequence"},
        "domain": domain, "min_coverage": min_coverage, "shared_pdb_ids": shared, "shared_pdb_count": len(shared),
        "domain_structure_count": len(domain_rows), "unmapped_shared_pdb_count": unmapped, "domain_structures": domain_rows,
        "selected_pdb": selected, "selected_pdb_metadata": protein["pdb"].get(selected) if selected else None,
        "selected_reason": ("highest coverage then best resolution among domain-overlapping PDB entries" if domain_rows and not pdb_id else ("user-supplied" if pdb_id else "alphabetical shared identifier; domain interval unavailable")),
        "caveat": "PROSITE.AUX 3D entries are profile-wide. Identifier overlap with UniProt PDB cross-references is not a residue-level mapping. Domain structures are those whose UniProt Chains interval overlaps the M-pro coordinates on this accession. That overlap is still an annotation agreement, not experimental proof of activity, and it is not an independent ScanProsite match.",
        "sources": {"prosite_primary": f"{BASES['prosite']}/{PROSITE_ID}.txt", "prosite_aux": f"{BASES['prosite']}/{PROSITE_ID}.aux", "uniprot": f"{BASES['uniprot']}/{accession}.json"},
    }
    if selected and include_structure:
        url = f"{BASES['pdb']}/{selected}.cif"
        report["structure"] = parse_cif(download(session, url), selected)
        report["sources"]["cif"] = url
    fasta = SeqRecord(Seq(protein["sequence"]), id=accession,
                      description=f"{protein['protein_name'] or 'protein'} | {protein['organism'] or 'unknown organism'}")
    return report, fasta, domain_structures_csv(domain_rows)
