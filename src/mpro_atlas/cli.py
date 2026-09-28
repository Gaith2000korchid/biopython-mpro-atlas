import argparse
import json
from pathlib import Path
from Bio import SeqIO
from .core import DEFAULT_ACCESSION, build_report, session_with_retries


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trace coronavirus M-pro PROSITE PS51442 to UniProt and PDB, then keep structures whose annotated chains overlap the M-pro interval"
    )
    parser.add_argument("--accession", default=DEFAULT_ACCESSION, help="Exact UniProtKB accession")
    parser.add_argument("--pdb", help="Choose a PDB cross-reference that overlaps the M-pro interval")
    parser.add_argument("--domain-start", type=int, help="Inclusive UniProt start of the M-pro interval")
    parser.add_argument("--domain-end", type=int, help="Inclusive UniProt end of the M-pro interval")
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=0.0,
        help="Minimum overlap with the M-pro interval, as a percentage of domain length (default: 0)",
    )
    parser.add_argument("--skip-structure", action="store_true", help="Skip potentially large mmCIF download")
    parser.add_argument("--out", type=Path, default=Path("results"))
    args = parser.parse_args()
    report, fasta, table = build_report(
        session_with_retries(),
        args.accession,
        args.pdb,
        not args.skip_structure,
        args.domain_start,
        args.domain_end,
        args.min_coverage,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (args.out / "domain_structures.csv").write_text(table, encoding="utf-8")
    with (args.out / "protein.fasta").open("w", encoding="utf-8") as handle:
        SeqIO.write(fasta, handle, "fasta")
    print(
        f"Saved {args.out / 'report.json'}, {args.out / 'domain_structures.csv'} "
        f"and {args.out / 'protein.fasta'} "
        f"({report['domain_structure_count']} domain structures / {report['shared_pdb_count']} shared IDs)"
    )


if __name__ == "__main__":
    main()
