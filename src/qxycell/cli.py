"""Command-line wrapper around the Python API."""

from __future__ import annotations

import argparse

from qxycell.checks import check
from qxycell.pipeline import import_cells


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qxycell")
    subparsers = parser.add_subparsers(dest="command", required=True)

    check_parser = subparsers.add_parser("check", help="Inspect and validate a QuPath project folder.")
    check_parser.add_argument("project_dir", help="Path to the QuPath project folder.")
    check_parser.add_argument(
        "--out",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Output folder for check reports. Defaults to <project-name>_check_YYMMDD_HHMM beside project_dir.",
    )
    check_parser.add_argument(
        "--count-rows",
        action="store_true",
        help="Count measurement table rows during check. Slower for large tables.",
    )
    check_parser.add_argument(
        "--threshold-file",
        default=None,
        help="Explicit threshold TSV/CSV to validate instead of auto-selecting from the project folder.",
    )

    import_parser = subparsers.add_parser(
        "import-cells",
        help="Import QuPath cell measurements into AnnData.",
    )
    import_parser.add_argument("project_dir", help="Path to the QuPath project folder.")
    import_parser.add_argument(
        "--out",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Output folder for AnnData and tables. Defaults to <project-name>_run_YYMMDD_HHMM beside project_dir.",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.command == "check":
        report = check(
            args.project_dir,
            output_dir=args.output_dir,
            count_rows=args.count_rows,
            threshold_file=args.threshold_file,
        )
        status = "PASS" if report.ok else "FAIL"
        print(f"QXYCell check {status}")
        print(f"Report: {report.output_dir / 'check_report.txt'}")
        print(f"Measurement files: {len(report.measurement_files)}")
        print(f"Classifier definitions: {len(report.classifiers)}")
        print(f"GeoJSON files: {len(report.geojson_files)}")
        print(f"Errors: {report.n_errors}; warnings: {report.n_warnings}")
        return 0 if report.ok else 1

    if args.command == "import-cells":
        adata = import_cells(
            args.project_dir,
            output_dir=args.output_dir,
        )
        print("QXYCell cell import complete")
        print(f"H5AD: {adata.uns['qxycell']['h5ad_path']}")
        print(f"Output: {adata.uns['qxycell']['run_dir']}")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
