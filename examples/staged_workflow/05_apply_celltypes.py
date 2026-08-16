"""Stage 5: replace cell-type outputs using the reviewed YAML logic."""

import qxycell as qxy

from config import CELLTYPE_YAML, OUTPUT_DIR


def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.celltype(adata, CELLTYPE_YAML)


if __name__ == "__main__":
    main()
