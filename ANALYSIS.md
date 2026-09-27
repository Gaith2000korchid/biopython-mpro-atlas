# Analysis: what did the empty PDB list mean?

## Question and method

The Coursera notebook called `ExPASy.get_prosite_raw('PS51442')`, parsed it with Biopython and obtained `record.pdb_structs == []`. The hypothesis that all PDB associations had disappeared can be tested against the PROSITE format documentation and the entry's auxiliary file.

We parsed `PS51442.txt` with `Bio.ExPASy.Prosite`, extracted the `3D` and `DR` lines of `PS51442.aux`, fetched **one** UniProtKB entry (`P0DTD1`), compared the accession's PDB cross-references with the profile's PDB list, then retrieved `6LU7.cif` for an example structure. The report keeps source URLs and retrieval time.

## Observation from 27 September 2026 UTC

| Evidence | Result | Scope |
| --- | ---: | --- |
| `PS51442.txt` / `record.pdb_structs` | 0 | Main profile record only |
| `PS51442.aux` / `3D` identifiers | 2,089 | Entire profile, many coronavirus proteins |
| UniProtKB `P0DTD1` / PDB cross-references | 3,363 | Entire SARS-CoV-2 replicase polyprotein 1ab |
| IDs shared by the two lists | 1,577 | Identifier overlap, not verified chain coverage |
| `6LU7` UniProt chain annotation | A=3264–3569 | Entry's position range |
| `6LU7` mmCIF parsed by Biopython | Chain A: 306 modeled residues; chain C: 3 | Observed residue counts |

The range 3264–3569 has 306 positions, consistent with chain A's residue count. These counts are time sensitive, and the accession covers an entire polyprotein: **a UniProt PDB cross-reference alone is not a claim that every PDB entry covers M-pro**. A residue-level mapping should precede any large-scale structural analysis. The 2,089 profile entries must not be mislabeled "the PDB structures of P0DTD1".

## Interpretation

The [PROSITE user manual](https://prosite.expasy.org/prosuser.html) documents that `3D` lines have been located in `PROSITE.AUX` since release 2024_03. [PS51442](https://prosite.expasy.org/PS51442) is the coronavirus M-pro domain profile, and its [auxiliary file](https://prosite.expasy.org/PS51442.aux) still contains many `3D` lines. This is a change in which file stores the annotation; the empty Biopython attribute describes the *primary record*. It is not evidence that PROSITE removed every PDB link.

The original notebook's broad UniProt `/search` request also yielded a list without checking subsequent pages. UniProt documents a default page size of 25 and a maximum size of 500. This project uses an exact accession instead, so its analysis has a defined biological unit.

## Next biological questions

Map individual PDB polymer entities to the M-pro region and inspect sequence identity, coverage, method, resolution and ligands. A genuine ScanProsite result, BLAST comparison, PubMed literature selection and KEGG annotation can then be added as **distinct evidence**, each with its own input, uncertainty and provenance.
