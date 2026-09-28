# Changelog

## 0.3.0 — 2026-09-28

Stable snapshot of the PS51442 / P0DTD1 case study.

- Default `--min-coverage` is 50. Coverage tiers `any_overlap`, `min_50` and `min_95` are reported.
- Automatic selection no longer falls back to `shared[0]`. If nothing meets the threshold, `selected_pdb` is empty.
- PDBe SIFTS UniProt mappings are applied to the selected structure (`--skip-sifts`, `--sifts-limit`).
- Report provenance includes SHA-256 hashes of the raw PROSITE and UniProt payloads.
- Offline tests cover AUX parsing, domain filtering, SIFTS missing-residue math and the missing fallback.

Live snapshot on 28 September 2026 UTC: 2,089 AUX PDB IDs, 1,577 shared with P0DTD1, 1,575 at ≥50% of 3264–3569. Rejected examples: `8H7K`, `9EZ6`.

## 0.2.0 — 2026-09-28

- Parse UniProt `Chains` ranges and keep structures that overlap the M-pro interval.
- Rank remaining entries by coverage then resolution.
- Write `domain_structures.csv`.

## 0.1.0 — 2026-09-27

- Explain empty `pdb_structs` on the primary PROSITE record.
- Parse `PS51442.aux` and intersect identifiers with one UniProt accession.
