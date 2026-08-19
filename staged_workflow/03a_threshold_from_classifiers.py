"""Stage 3A: replace marker positivity using classifier JSON thresholds only."""

import qxycell as qxy

from config import OUTPUT_DIR, PROJECT_DIR


def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.threshold_from_classifiers(adata, project_dir=PROJECT_DIR)


if __name__ == "__main__":
    main()

