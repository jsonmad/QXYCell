"""Stage 4: generate the LLM prompt used to draft cell-type YAML."""

import qxycell as qxy

from config import CELLTYPE_CONTEXT, OUTPUT_DIR


def main() -> None:
    adata = qxy.load(OUTPUT_DIR)
    qxy.celltype_prompt(adata, context=CELLTYPE_CONTEXT)


if __name__ == "__main__":
    main()

