# Recorded source evidence

`report-summary.json` is the historical 28 September 2026 summary. Its original source payloads were not retained here; new fingerprints must not be attributed retroactively to that run.

`report-summary-2026-10-09.json` records a new live execution with v0.3.1, including the source manifest and SHA-256 fingerprints. `sources-2026-10-09.zip` contains the five successfully retrieved decoded UTF-8 payloads and their manifest: primary PROSITE, PROSITE.AUX, UniProt, PDBe/SIFTS and the selected RCSB mmCIF. These are the texts consumed by the parsers, not original compressed HTTP wire bytes.

## Replay the dated case without network access

After installing the frozen dependencies and the local package, from the repository root:

```bash
python -m zipfile -e examples/sources-2026-10-09.zip results/snapshot-2026-10-09
mpro-atlas --accession P0DTD1 --pdb 6LU7 --offline-sources results/snapshot-2026-10-09/sources --out results/replayed-2026-10-09
```

Offline mode validates payload checksums, refuses missing URLs and identifies itself as `offline_snapshot` in the report. Its replay timestamp is different from the original retrieval date, which remains in the manifest. Biological counts and mappings can be compared using the same recorded inputs and locked parser dependencies.

To record a new live execution, use a new output directory:

```bash
mpro-atlas --accession P0DTD1 --pdb 6LU7 --out results/new-live-run
```

A live run refuses to overwrite a directory containing a source manifest. Failed HTTP responses are not archived as successful evidence. A partial source snapshot can exist if a later request fails; do not call it a completed analysis.

## Attribution

Source records originate from SIB PROSITE, UniProtKB, PDBe/SIFTS and RCSB PDB. The manifest and report retain their URLs. Their data/source terms apply independently of this repository's MIT code license. No Coursera source material or private health data is included.

Annotations and segment mappings do not demonstrate enzymatic activity or per-atom occupancy. The archived data freezes a case study, not the current state of the source databases.
