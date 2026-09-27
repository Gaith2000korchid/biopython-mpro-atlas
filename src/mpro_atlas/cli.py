import argparse
import json
from pathlib import Path
from Bio import SeqIO
from .core import DEFAULT_ACCESSION, build_report, session_with_retries


def main() -> None:
    parser = argparse.ArgumentParser(description="Trace coronavirus M-pro PROSITE PS51442 to UniProt and PDB")
    parser.add_argument("--accession", default=DEFAULT_ACCESSION, help="Exact UniProtKB accession")
    parser.add_argument("--pdb", help="Choose a PDB cross-reference from the UniProt entry")
    parser.add_argument("--skip-structure", action="store_true", help="Skip potentially large mmCIF download")
    parser.add_argument("--out", type=Path, default=Path("results"))
    args = parser.parse_args()
    report, fasta = build_report(session_with_retries(), args.accession, args.pdb, not args.skip_structure)
    args.out.mkdir(parents=True, exist_ok=True)
    (args.out / "report.json").write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (args.out / "protein.fasta").open("w", encoding="utf-8") as handle:
        SeqIO.write(fasta, handle, "fasta")
    print(f"Saved {args.out / 'report.json'} and {args.out / 'protein.fasta'}")


if __name__ == "__main__":
    main()
