# Coronavirus M-pro: from PROSITE profile to protein structures

An independently built, reproducible bioinformatics case study inspired by Coursera's **Access Bioinformatics Databases with Biopython** Guided Project. It investigates why `Bio.ExPASy.Prosite.read(...).pdb_structs` returns `[]` for `PS51442` even though experimental structures exist.

**Finding:** `PS51442` is the *coronavirus main protease (M-pro)* profile. Since PROSITE release 2024_03, `3D` structure references reside in the separate `PROSITE.AUX` entry (`PS51442.aux`), while the regular `PS51442.txt` record contains the profile definition. An empty `pdb_structs` on the regular record is therefore a record-format issue, not evidence that the protein has no PDB structures. See [PROSITE record](https://prosite.expasy.org/PS51442), [auxiliary record](https://prosite.expasy.org/PS51442.aux), and [PROSITE user manual](https://prosite.expasy.org/prosuser.html).

## Reproduce

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install -e '.[test]'
python -m pytest -q
mpro-atlas --accession P0DTD1 --pdb 6LU7 --out results
```

The command contacts the live PROSITE, UniProtKB and RCSB PDB services. It creates `results/report.json` and `results/protein.fasta`; generated outputs are excluded from git because external databases change. `--skip-structure` skips the mmCIF download; `--pdb 6LU7` requests a specific PDB ID **provided UniProt currently cross-references it**. Use `mpro-atlas --help` for options. Network access is needed only for the command, not for tests.

The checked-in [`examples/report-summary.json`](examples/report-summary.json) records an actual run on 27 September 2026 UTC: `pdb_structs=[]` in the primary record, 2,089 profile-wide IDs in `PS51442.aux`, and 6LU7 mapped by UniProt to polyprotein positions 3264–3569. Chain A of the retrieved structure contains 306 modeled residues. Counts depend on the release and will change. The example is deliberately compact; rerun the command to obtain the complete current cross-reference lists.

See [ANALYSIS.md](ANALYSIS.md) for the evidence table, biological interpretation and remaining validation work.

### Windows PowerShell

```powershell
py -m venv .venv
.\.venv\Scripts\python.exe -m pip install -e ".[test]"
.\.venv\Scripts\python.exe -m pytest -q
.\.venv\Scripts\mpro-atlas.exe --accession P0DTD1 --pdb 6LU7 --out results
```

## How the evidence is joined

1. Biopython parses the primary PROSITE profile `PS51442.txt`.
2. A small parser reads `PS51442.aux` for the *profile-wide* PDB list and the matched UniProtKB accessions.
3. The tool fetches **one exact UniProtKB accession** (default `P0DTD1`, SARS-CoV-2 replicase polyprotein 1ab), extracts its protein sequence, PROSITE references and PDB cross-references. This avoids a truncated search result and accidental mixing of proteins from different organisms.
4. It intersects PROSITE.AUX PDB IDs with that protein's PDB IDs and parses one experimental structure from mmCIF using `Bio.PDB.MMCIFParser`.

The intersection is an **identifier-level cross-reference**, not proof that every chain covers the M-pro region. UniProt's `Chains` metadata and residue mappings must be checked for any residue-level biological claim. This tool does not perform a new PROSITE sequence scan, BLAST search or KEGG pathway analysis; those are separate questions and are not implied by the PDB cross-references. `report.json` records the source URLs and retrieval time so results can be audited.

## Known limits

- Live endpoints may be unavailable, rate limited or updated; offline tests validate parsing and report logic, not remote service availability.
- The CLI selects the alphabetically first shared PDB ID unless `--pdb` is supplied; this is **not** a quality ranking. Inspect method, resolution, chain coverage and ligands before drawing structural conclusions.
- The accession `P0DTD1` describes a *polyprotein*, not an isolated mature M-pro sequence.
- For a broad UniProt search such as `xref:prosite-PS51442`, the `/search` endpoint is paginated (default 25, max 500 per page). A single response must not be reported as the total count. See [UniProt API paper](https://academic.oup.com/nar/article/53/W1/W547/8126256).

## Attribution

Concept inspired by the named Coursera Guided Project. All implementation, analysis and offline test fixtures here are newly authored; Coursera course files and videos are not redistributed. Data and annotations remain attributed to [SIB PROSITE](https://prosite.expasy.org/), [UniProtKB](https://www.uniprot.org/) and [RCSB PDB](https://www.rcsb.org/).
