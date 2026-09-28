import argparse
import json
from pathlib import Path
from Bio import SeqIO
from .core import DEFAULT_ACCESSION, DEFAULT_MIN_COVERAGE, build_report, session_with_retries


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Trace PS51442 to UniProt/PDB and keep structures whose chains cover the M-pro interval"
    )
    parser.add_argument("--accession", default=DEFAULT_ACCESSION, help="Exact UniProtKB accession")
    parser.add_argument("--pdb", help="Choose a PDB cross-reference that meets min-coverage")
    parser.add_argument("--domain-start", type=int, help="Inclusive UniProt start of the M-pro interval")
    parser.add_argument("--domain-end", type=int, help="Inclusive UniProt end of the M-pro interval")
    parser.add_argument(
        "--min-coverage",
        type=float,
        default=DEFAULT_MIN_COVERAGE,
        help="Minimum UniProt Chains overlap with the M-pro interval (default: 50)",
    )
    parser.add_argument("--skip-structure", action="store_true", help="Skip mmCIF download")
    parser.add_argument("--skip-sifts", action="store_true", help="Skip PDBe/SIFTS mapping")
    parser.add_argument(
        "--sifts-limit",
        type=int,
        default=1,
        help="Max number of PDB entries to enrich with SIFTS (default: 1, the selected structure)",
    )
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
        not args.skip_sifts,
        args.sifts_limit,
    )
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    (args.out / "domain_structures.csv").write_text(table, encoding="utf-8")
    with (args.out / "protein.fasta").open("w", encoding="utf-8") as handle:
        SeqIO.write(fasta, handle, "fasta")
    selected = report["selected_pdb"] or "none"
    tiers = report.get("coverage_tiers") or {}
    print(
        f"Saved {args.out / 'report.json'}, {args.out / 'domain_structures.csv'} "
        f"and {args.out / 'protein.fasta'} "
        f"(selected={selected}; >=50%={tiers.get('min_50')} / shared={report['shared_pdb_count']})"
    )


if __name__ == "__main__":
    main()
