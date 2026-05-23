"""Command-line wrapper around the Python API."""

from __future__ import annotations

import argparse
from pathlib import Path

from quxycell.checks import check
from quxycell.pipeline import run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="quxycell")
    subparsers = parser.add_subparsers(dest="command", required=True)

    check_parser = subparsers.add_parser("check", help="Inspect and validate a QuPath export.")
    check_parser.add_argument("project_dir", help="Path to a manually exported QuPath project folder.")
    check_parser.add_argument(
        "--out",
        "--output-dir",
        dest="output_dir",
        default="quxycell_output",
        help="Output folder for check reports.",
    )
    check_parser.add_argument(
        "--count-rows",
        action="store_true",
        help="Count measurement table rows during check. Slower for large exports.",
    )

    run_parser = subparsers.add_parser("run", help="Run QUXYCell on a QuPath export.")
    run_parser.add_argument("project_dir", help="Path to a manually exported QuPath project folder.")
    run_parser.add_argument(
        "--out",
        "--output-dir",
        dest="output_dir",
        default="quxycell_output",
        help="Output folder for QUXYCell reports and AnnData.",
    )
    run_parser.add_argument(
        "--pixel-size-um",
        type=float,
        default=0.28,
        help="Pixel size used to scale QuPath GeoJSON coordinates into microns.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "check":
        report = check(args.project_dir, output_dir=args.output_dir, count_rows=args.count_rows)
        status = "PASS" if report.ok else "FAIL"
        print(f"QUXYCell check {status}")
        print(f"Report: {Path(args.output_dir).expanduser().resolve() / 'check_report.txt'}")
        print(f"Measurement files: {len(report.measurement_files)}")
        print(f"Classifier JSON files: {len(report.classifiers)}")
        print(f"GeoJSON files: {len(report.geojson_files)}")
        print(f"Errors: {report.n_errors}; warnings: {report.n_warnings}")
        return 0 if report.ok else 1

    if args.command == "run":
        result = run(
            args.project_dir,
            output_dir=args.output_dir,
            pixel_size_um=args.pixel_size_um,
        )
        print("QUXYCell run complete")
        print(f"H5AD: {result.h5ad_path}")
        print(f"Output: {result.output_dir}")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
