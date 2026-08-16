"""Stage 1: create the base AnnData checkpoint from QuPath measurements."""

import qxycell as qxy

from config import OUTPUT_DIR, PROJECT_DIR


def main() -> None:
    qxy.import_measurements(PROJECT_DIR, output_dir=OUTPUT_DIR)


if __name__ == "__main__":
    main()

