# Coronavirus M-pro: from PROSITE profile to protein structures

An independently built, reproducible bioinformatics case study inspired by Coursera's **Access Bioinformatics Databases with Biopython** Guided Project. It investigates why `Bio.ExPASy.Prosite.read(...).pdb_structs` returns `[]` for `PS51442` even though experimental structures exist, then keeps only the PDB entries whose UniProt chain coordinates cover at least 50% of the M-pro interval, and optionally checks the selected structure with PDBe/SIFTS.

**Finding:** `PS51442` is the *coronavirus main protease (M-pro)* profile. Since PROSITE release 2024_03, `3D` structure references reside in the separate `PROSITE.AUX` entry (`PS51442.aux`), while the regular `PS51442.txt` record contains the profile definition. An empty `pdb_structs` on the regular record is therefore a record-format issue, not evidence that the protein has no PDB structures. See [PROSITE record](https://prosite.expasy.org/PS51442), [auxiliary record](https://prosite.expasy.org/PS51442.aux), and [PROSITE user manual](https://prosite.expasy.org/prosuser.html).

Identifier overlap is still not enough. `P0DTD1` is the SARS-CoV-2 replicase polyprotein 1ab. A PDB cross-reference on that accession may cover nsp5 / M-pro or another region. Version 0.2 parses UniProt `Chains` annotations and compares them with the M-pro interval. Version 0.3 raises the default coverage threshold to 50%, refuses to fall back to an unfiltered shared PDB ID, and maps the selected structure through PDBe SIFTS segments.

## Reproduce

```bash
python -m venv .venv
source .venv/bin/activate  # Windows PowerShell: .\.venv\Scripts\Activate.ps1
python -m pip install -e '.[test]'
python -m pytest -q
mpro-atlas --accession P0DTD1 --pdb 6LU7 --out results
```

The command contacts the live PROSITE, UniProtKB, PDBe and RCSB PDB services. It creates `results/report.json`, `results/domain_structures.csv` and `results/protein.fasta`; generated outputs are excluded from git because external databases change. `--skip-structure` skips the mmCIF download; `--skip-sifts` skips PDBe mapping; `--sifts-limit` controls how many passing PDB IDs are enriched. `--pdb 6LU7` requests a specific PDB ID **provided UniProt currently cross-references it and the annotated chains meet `--min-coverage` on the M-pro interval**. The default `--min-coverage` is 50. If no structure meets that threshold, the report leaves `selected_pdb` empty instead of picking `shared[0]`. Use `--domain-start` / `--domain-end` to override the interval. Network access is needed only for the command, not for tests.

The checked-in [`examples/report-summary.json`](examples/report-summary.json) records a run on 28 September 2026 UTC: `pdb_structs=[]` in the primary record, 2,089 profile-wide IDs in `PS51442.aux`, 1,577 IDs shared with `P0DTD1`, and 1,575 of those covering at least 50% of 3264–3569. Two shared IDs fail that filter (`8H7K`, `9EZ6`). `6LU7` remains a valid example: UniProt maps chain A to 3264–3569, SIFTS maps the same UniProt interval to PDB residues 1–306, and the mmCIF contains 306 modeled residues on chain A. Counts depend on the release and will change.

See [ANALYSIS.md](ANALYSIS.md) for the evidence table, biological interpretation and remaining validation work.

### Windows PowerShell

```powershell
py -m venv .venv
.\ .venv\Scripts\python.exe -m pip install -e ".[test]"
