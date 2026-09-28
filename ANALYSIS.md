# Analysis: what did the empty PDB list mean?

## Question and method

The Coursera notebook called `ExPASy.get_prosite_raw('PS51442')`, parsed it with Biopython and obtained `record.pdb_structs == []`. The hypothesis that all PDB associations had disappeared can be tested against the PROSITE format documentation and the entry's auxiliary file.

We parsed `PS51442.txt` with `Bio.ExPASy.Prosite`, extracted the `3D` and `DR` lines of `PS51442.aux`, fetched **one** UniProtKB entry (`P0DTD1`), located the nsp5 / M-pro interval on that polyprotein, compared each shared PDB `Chains` annotation with that interval at a 50% coverage threshold, mapped the selected structure through PDBe SIFTS, then retrieved `6LU7.cif`. The report keeps source URLs, retrieval time and payload hashes.

## Observation from 28 September 2026 UTC

| Evidence | Result | Scope |
| --- | ---: | --- |
| `PS51442.txt` / `record.pdb_structs` | 0 | Main profile record only |
| `PS51442.aux` / `3D` identifiers | 2,089 | Entire profile, many coronavirus proteins |
| UniProtKB `P0DTD1` / PDB cross-references | 3,363 | Entire SARS-CoV-2 replicase polyprotein 1ab |
| IDs shared by the two lists | 1,577 | Identifier overlap, not verified chain coverage |
| UniProt feature used as M-pro interval | 3264–3569 | Chain "3C-like proteinase nsp5" |
| Shared IDs whose `Chains` overlap that interval | 1,576 | Any positive overlap |
| Shared IDs with ≥50% domain coverage | 1,575 | Default annotation-level filter |
| Shared IDs rejected at 50% coverage | 2 | `8H7K` (B=3258–3267, 1.3%); `9EZ6` (C=6447–6457, 0%) |
| Default selection rule | no `shared[0]` fallback | Empty `selected_pdb` if nothing meets 50% |
| `6LU7` UniProt chain annotation | A=3264–3569 | Entry's position range |
| `6LU7` PDBe SIFTS UniProt segment | A unp 3264–3569 → PDB 1–306 | Segment alignment |
| `6LU7` mmCIF parsed by Biopython | Chain A: 306 modeled residues; chain C: 3 | Observed residue counts |

The range 3264–3569 has 306 positions, consistent with chain A's residue count. These counts are time sensitive. The 2,089 profile entries must not be mislabeled "the PDB structures of P0DTD1". The 1,577 shared IDs must not be reported as 1,577 independent proofs of M-pro coverage: two of them fail a simple range test.

On this accession and this PROSITE profile, identifier intersection already removes most non-nsp5 structures of `P0DTD1`. The residue-range step is still required. It converts an implicit assumption into a stated predicate, ranks remaining entries by coverage and resolution, and exposes the exceptions.

## Interpretation

The [PROSITE user manual](https://prosite.expasy.org/prosuser.html) documents that `3D` lines have been located in `PROSITE.AUX` since release 2024_03. [PS51442](https://prosite.expasy.org/PS51442) is the coronavirus M-pro domain profile, and its [auxiliary file](https://prosite.expasy.org/PS51442.aux) still contains many `3D` lines. This is a change in which file stores the annotation; the empty Biopython attribute describes the *primary record*. It is not evidence that PROSITE removed every PDB link.

The original notebook's broad UniProt `/search` request also yielded a list without checking subsequent pages. UniProt documents a default page size of 25 and a maximum size of 500. This project uses an exact accession instead, so its analysis has a defined biological unit.

UniProt `Chains` coordinates remain annotations. Version 0.3 adds a PDBe SIFTS check on the selected structure: for `6LU7` the UniProt interval 3264–3569 is aligned to PDB residues 1–306 with no missing domain residues in the merged segments. That is still a segment alignment, not a new experimental result and not a per-atom occupancy proof.

## Remaining biological questions

SIFTS is now applied to the selected structure (and optionally a short tail via `--sifts-limit`). A residue-occupancy audit of missing atoms inside those segments is still out of scope. Inspect ligands and mutations on the filtered table. A genuine ScanProsite result, BLAST comparison, PubMed literature selection and KEGG annotation can then be added as **distinct evidence**, each with its own input, uncertainty and provenance.
