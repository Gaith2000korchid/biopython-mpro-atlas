# Coronavirus M-pro: from PROSITE profile to protein structures

An independently built, reproducible bioinformatics case study inspired by Coursera's **Access Bioinformatics Databases with Biopython** Guided Project. It investigates why `Bio.ExPASy.Prosite.read(...).pdb_structs` returns `[]` for `PS51442` even though experimental structures exist, then keeps only the PDB entries whose UniProt chain coordinates overlap the M-pro interval.

**Finding:** `PS51442` is the *coronavirus main protease (M-pro)* profile. Since PROSITE release 2024_03, `3D` structure references reside in the separate `PROSITE.AUX` entry (`PS51442.aux`), while the regular `PS51442.txt` record contains the profile definition. An empty `pdb_structs` on the regular record is therefore a record-format issue, not evidence that the protein has no PDB structures. See [PROSITE record](https://prosite.expasy.org/PS51442), [auxiliary record](https://prosite.expasy.org/PS51442.aux), and [PROSITE user manual](https://prosite.expasy.org/prosuser.html).

Identifier overlap is still not enough. `P0DTD1` is the SARS-CoV-2 replicase polyprotein 1ab. A PDB cross-reference on that accession may cover nsp5 / M-pro or another region. Version 0.2 parses UniProt `Chains` annotations, compares them with the M-pro interval, and writes a table of domain-overlapping structures.

## Reproduce

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install -e '.[test]'
python -m pytest -q
mpro-atlas --accession P0DTD1 --pdb 6LU7 --out results
```

The command contacts the live PROSITE, UniProtKB and RCSB PDB services. It creates `results/report.json`, `results/domain_structures.csv` and `results/protein.fasta`; generated outputs are excluded from git because external databases change. `--skip-structure` skips the mmCIF download; `--pdb 6LU7` requests a specific PDB ID **provided UniProt currently cross-references it and the annotated chains overlap the M-pro interval**. Use `--domain-start` / `--domain-end` to override the interval and `--min-coverage` to require a minimum fraction of the domain. Use `mpro-atlas --help` for options. Network access is needed only for the command, not for tests.

The checked-in [`examples/report-summary.json`](examples/report-summary.json) records an actual run on 28 September 2026 UTC: `pdb_structs=[]` in the primary record, 2,089 profile-wide IDs in `PS51442.aux`, 1,577 IDs shared with `P0DTD1`, and 1,575 of those overlapping at least 50% of the nsp5 interval 3264–3569. Two shared IDs fail that filter (`8H7K`, `9EZ6`). `6LU7` remains a valid example: UniProt maps chain A to 3264–3569 and the mmCIF contains 306 modeled residues on chain A. Counts depend on the release and will change. The example is deliberately compact; rerun the command to obtain the complete current table.

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
2. A small parser reads `PS51442.aux` for the *profile-wide* PDB list and the matched UniProtKB accessions. Tokens that are not PDB IDs are ignored.
3. The tool fetches **one exact UniProtKB accession** (default `P0DTD1`, SARS-CoV-2 replicase polyprotein 1ab), extracts its protein sequence, domain/chain features, PROSITE references and PDB cross-references. This avoids a truncated search result and accidental mixing of proteins from different organisms.
4. It locates the M-pro interval on that accession (UniProt feature linked to `PS51442`, otherwise a description match such as "3C-like proteinase nsp5", otherwise the documented P0DTD1 default 3264–3569).
5. It intersects PROSITE.AUX PDB IDs with that protein's PDB IDs, then **keeps only entries whose UniProt `Chains` range overlaps the M-pro interval**.
6. It ranks the remaining structures by coverage, then by resolution, and can parse one experimental structure from mmCIF using `Bio.PDB.MMCIFParser`.

The first intersection is an **identifier-level cross-reference**. The second step is an **annotation-level residue range comparison** on the same UniProt entry. Neither is proof of enzymatic activity, and the tool does not perform a new PROSITE sequence scan, BLAST search or KEGG pathway analysis. `report.json` records the source URLs and retrieval time so results can be audited.

## Known limits

- Live endpoints may be unavailable, rate limited or updated; offline tests validate parsing and report logic, not remote service availability.
- UniProt `Chains` coordinates are used as-is. They are not a SIFTS residue-by-residue alignment and can disagree with author numbering in the mmCIF.
- Without `--pdb`, the CLI selects the domain-overlapping entry with the highest coverage and, at equal coverage, the best reported resolution. Inspect method, ligands and mutations before drawing structural conclusions.
- The accession `P0DTD1` describes a *polyprotein*, not an isolated mature M-pro sequence. The interval 3264–3569 is the mature nsp5 region on that polyprotein.
- For a broad UniProt search such as `xref:prosite-PS51442`, the `/search` endpoint is paginated (default 25, max 500 per page). A single response must not be reported as the total count. See [UniProt API paper](https://academic.oup.com/nar/article/53/W1/W547/8126256).

## Attribution

Concept inspired by the named Coursera Guided Project. All implementation, analysis and offline test fixtures here are newly authored; Coursera course files and videos are not redistributed. Data and annotations remain attributed to [SIB PROSITE](https://prosite.expasy.org/), [UniProtKB](https://www.uniprot.org/) and [RCSB PDB](https://www.rcsb.org/).
