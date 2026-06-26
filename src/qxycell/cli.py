"""Command-line wrapper around the Python API."""

from __future__ import annotations

import argparse

from qxycell.checks import check
from qxycell.pipeline import run


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="qxycell")
    subparsers = parser.add_subparsers(dest="command", required=True)

    check_parser = subparsers.add_parser("check", help="Inspect and validate a QuPath export.")
    check_parser.add_argument("project_dir", help="Path to a manually exported QuPath project folder.")
    check_parser.add_argument(
        "--out",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Output folder for check reports. Defaults to qxy_outputs_YYMMDD-HHMM.",
    )
    check_parser.add_argument(
        "--count-rows",
        action="store_true",
        help="Count measurement table rows during check. Slower for large exports.",
    )
    check_parser.add_argument(
        "--threshold-file",
        default=None,
        help="Explicit threshold TSV/CSV to validate instead of auto-selecting from the project folder.",
    )

    run_parser = subparsers.add_parser("run", help="Run QXYCell on a QuPath export.")
    run_parser.add_argument("project_dir", help="Path to a manually exported QuPath project folder.")
    run_parser.add_argument(
        "--out",
        "--output-dir",
        dest="output_dir",
        default=None,
        help="Output folder for QXYCell reports and AnnData. Defaults to qxy_outputs_YYMMDD-HHMM.",
    )
    run_parser.add_argument(
        "--pixel-size-um",
        type=float,
        default=0.28,
        help="Pixel size used to scale QuPath GeoJSON coordinates into microns.",
    )
    run_parser.add_argument(
        "--celltype-logic",
        default=None,
        help="Optional YAML file with ordered cell type rules. Implies --apply-thresholds.",
    )
    run_parser.add_argument(
        "--threshold-file",
        default=None,
        help="Explicit threshold TSV/CSV to use when --apply-thresholds is set.",
    )
    run_parser.add_argument(
        "--apply-thresholds",
        action="store_true",
        help="Apply threshold definitions after importing AnnData to create <marker>_pos columns.",
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

    if args.command == "run":
        adata = run(
            args.project_dir,
            output_dir=args.output_dir,
            pixel_size_um=args.pixel_size_um,
            celltype_logic=args.celltype_logic,
            threshold_file=args.threshold_file,
            apply_thresholds=args.apply_thresholds or args.celltype_logic is not None,
        )
        print("QXYCell run complete")
        print(f"H5AD: {adata.uns['qxycell']['h5ad_path']}")
        print(f"Output: {adata.uns['qxycell']['run_dir']}")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


if __name__ == "__main__":
    raise SystemExit(main())
