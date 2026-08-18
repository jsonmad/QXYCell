"""Stage 2b: remove cells inside annotations whose names contain configured text."""

import qxycell as qxy

from config import IGNORE_ANNOTATION_TEXT, OUTPUT_DIR


def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.remove_cells(
        adata,
        annotation_prefix="annotation__",
        remove_cells=IGNORE_ANNOTATION_TEXT,
        copy=False,
        verbose=True,
    )
    qxy.save(adata, output_dir=OUTPUT_DIR, verbose=True)

    tables_dir = OUTPUT_DIR / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)
    adata.obs.to_csv(tables_dir / "cells_obs.csv")


if __name__ == "__main__":
    main()
