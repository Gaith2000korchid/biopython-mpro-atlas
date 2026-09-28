"""Deterministic domain and SIFTS mapping helpers."""
from __future__ import annotations
from typing import Any
import json
import re
import requests

_CHAIN_ASSIGNMENT = re.compile(r"([A-Za-z0-9][A-Za-z0-9/]*)\s*=\s*(\d+)\s*-\s*(\d+)")


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


def filter_domain_structures(pdb_annotations, domain_start, domain_end, candidate_ids=None, min_coverage=50.0):
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


def coverage_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    return {
        "any_overlap": sum(1 for row in rows if row["coverage_pct"] > 0),
        "min_50": sum(1 for row in rows if row["coverage_pct"] >= 50),
        "min_95": sum(1 for row in rows if row["coverage_pct"] >= 95),
    }


def parse_sifts_uniprot(payload: dict[str, Any], pdb_id: str, accession: str) -> list[dict[str, Any]]:
    block = payload.get(pdb_id.lower()) or payload.get(pdb_id) or payload.get(pdb_id.upper()) or {}
    entry = ((block.get("UniProt") or {}).get(accession) or {})
    segments = []
    for mapping in entry.get("mappings") or []:
        start, end = mapping.get("start") or {}, mapping.get("end") or {}
        unp_start, unp_end = mapping.get("unp_start"), mapping.get("unp_end")
        if unp_start is None or unp_end is None:
            continue
        unp_start, unp_end = int(unp_start), int(unp_end)
        if unp_start > unp_end:
            unp_start, unp_end = unp_end, unp_start
        segments.append({
            "chain_id": mapping.get("chain_id"), "unp_start": unp_start, "unp_end": unp_end,
            "pdb_start": start.get("residue_number"), "pdb_end": end.get("residue_number"),
            "author_start": start.get("author_residue_number"), "author_end": end.get("author_residue_number"),
            "identity": mapping.get("identity"),
        })
    return segments


def sifts_domain_coverage(segments: list[dict[str, Any]], domain_start: int, domain_end: int) -> dict[str, Any]:
    intervals = []
    for segment in segments:
        lo = max(segment["unp_start"], domain_start)
        hi = min(segment["unp_end"], domain_end)
        if lo <= hi:
            intervals.append([lo, hi])
    intervals.sort()
    merged = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1] + 1:
            merged.append([lo, hi])
        else:
            merged[-1][1] = max(merged[-1][1], hi)
    covered = sum(hi - lo + 1 for lo, hi in merged)
    domain_length = domain_end - domain_start + 1
    missing, cursor = [], domain_start
    for lo, hi in merged:
        if cursor < lo:
            missing.append({"start": cursor, "end": lo - 1})
        cursor = hi + 1
    if cursor <= domain_end:
        missing.append({"start": cursor, "end": domain_end})
    return {
        "segments": segments,
        "overlap_residues": covered,
        "coverage_pct": round(100.0 * covered / domain_length, 2) if domain_length else 0.0,
        "missing_residues": sum(item["end"] - item["start"] + 1 for item in missing),
        "missing_ranges": missing,
        "mapping_source": "pdbe_sifts_uniprot",
    }


def fetch_sifts(session: requests.Session, pdb_id: str, pdbe_base: str) -> tuple[dict[str, Any], str]:
    url = f"{pdbe_base}/mappings/uniprot/{pdb_id.lower()}"
    response = session.get(url, timeout=30)
    response.raise_for_status()
    return json.loads(response.text), url


def attach_sifts(row: dict[str, Any], payload: dict[str, Any], accession: str, domain_start: int, domain_end: int) -> dict[str, Any]:
    segments = parse_sifts_uniprot(payload, row["pdb_id"], accession)
    row = dict(row)
    row["sifts"] = sifts_domain_coverage(segments, domain_start, domain_end) if segments else {
        "segments": [], "overlap_residues": 0, "coverage_pct": 0.0,
        "missing_residues": domain_end - domain_start + 1,
        "missing_ranges": [{"start": domain_start, "end": domain_end}],
        "mapping_source": "pdbe_sifts_uniprot",
        "note": "No SIFTS UniProt segment for this accession",
    }
    return row


def select_structure(pdb_id, protein, domain, domain_rows, min_coverage):
    if pdb_id:
        selected = pdb_id.upper()
        if selected not in protein["pdb"]:
            raise ValueError(f"{selected} is not cross-referenced by UniProt {protein['accession']}")
        if domain is None:
            return selected, "user-supplied; domain interval unavailable"
        row = map_pdb_to_domain(selected, protein["pdb"][selected], domain["start"], domain["end"])
        if row["coverage_pct"] < min_coverage:
            raise ValueError(
                f"{selected} does not meet min_coverage={min_coverage} on M-pro {domain['start']}-{domain['end']}"
            )
        return selected, "user-supplied"
    if domain_rows:
        return domain_rows[0]["pdb_id"], "highest UniProt-range coverage then best resolution among structures meeting min_coverage"
    if domain is not None:
        return None, "no PDB structure satisfies the requested domain coverage"
    return None, "domain interval unavailable; no structure selected"
