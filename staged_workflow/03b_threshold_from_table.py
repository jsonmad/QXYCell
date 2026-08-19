"""Stage 3B: replace marker positivity using one threshold table only."""

import qxycell as qxy

from config import OUTPUT_DIR, PROJECT_DIR, THRESHOLD_TABLE


def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.threshold_from_table(
        adata,
        THRESHOLD_TABLE,
        project_dir=PROJECT_DIR,
    )


if __name__ == "__main__":
    main()

