"""Stage 2: add or replace GeoJSON-derived annotations and cell polygons."""

import qxycell as qxy

from config import OUTPUT_DIR, PIXEL_SIZE_UM, PROJECT_DIR


def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.add_annotations(
        adata,
        project_dir=PROJECT_DIR,
        pixel_size_um=PIXEL_SIZE_UM,
    )


if __name__ == "__main__":
    main()

