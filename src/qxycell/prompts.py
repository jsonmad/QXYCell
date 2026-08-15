"""Prompt helpers for external LLM-assisted configuration drafting."""

from __future__ import annotations

from pathlib import Path
from textwrap import dedent

from qxycell.paths import resolve_output_dir
from qxycell.stage_state import complete_stage, prepare_stage


def _available_markers_from_adata(adata) -> list[str]:
    markers = []
    obs = getattr(adata, "obs", None)
    if obs is not None:
        for column in obs.columns:
            column = str(column)
            if column.endswith("_pos"):
                markers.append(column[:-4])

    return sorted(dict.fromkeys(marker for marker in markers if marker))


def celltype_prompt(
    adata,
    *,
    context: str | None = None,
    output_path: str | Path | None = None,
    save: bool = True,
    print_prompt: bool = True,
) -> str:
    """Create a copy-pasteable LLM prompt for drafting cell type logic YAML.

    The prompt asks an external LLM to create a first-pass QXYCell
    ``celltype_logic.yaml`` using only markers available in the loaded AnnData.
    """

    markers = _available_markers_from_adata(adata)
    if not markers:
        raise ValueError(
            "No thresholded markers were found in adata.obs '*_pos' columns. "
            "Run qxy.threshold() before creating a cell type prompt."
        )

    marker_lines = "\n".join(f"- {marker}" for marker in markers)
    context_text = context.strip() if context and context.strip() else "[No additional context provided]"

    prompt = (
        dedent(
            """
            You are helping draft a first-pass QXYCell cell type logic YAML file for spatial single-cell data imported from QuPath.

            QXYCell has already loaded the project into an AnnData object. Marker positivity columns exist in `adata.obs` as `<MARKER>_pos`. Use only the marker names provided below exactly as written.

            Available markers:
            """
        ).strip()
        + "\n"
        + marker_lines
        + "\n\n"
        + "Additional biological/project context:\n"
        + context_text
        + "\n\n"
        + dedent(
            """
            Return ONLY a valid YAML document inside a single fenced ```yaml code block.

            Critical formatting requirements:
            - Do not include explanatory text before or after the YAML.
            - Preserve YAML indentation exactly.
            - Do not convert YAML syntax into markdown bullet formatting.
            - YAML list items using `-` must remain literal plain-text YAML.
            - Do not use rich text formatting outside the YAML code block.

            Required structure:

            rules:
              - name: Example_celltype
                positive: [MARKER_A, MARKER_B]
                negative: [MARKER_C]

            features:
              Example_feature:
                positive: [MARKER_A]

            derived_features:
              Example_derived_feature:
                all_of: [Example_feature]

            Requirements:
            - Use only markers from the provided marker list.
            - Rules must be ordered from most specific to most general.
            - Rules are exclusive and first-match: the first matching rule assigns `celltype`.
            - Put rare/specific phenotypes before broad parent populations.
            - Use `positive` for required marker-positive calls.
            - Use `negative` for biologically useful exclusions and to prevent broad rules from capturing specific subtypes.
            - Use `any_positive` only when a cell type can be defined by any one of several markers.
            - Include cautious biologically sensible immune, stromal, endothelial, tumor, proliferation, and exhaustion states only if supported by the panel.
            - Do not invent markers, metadata fields, annotations, sample groups, treatment groups, or measurements.
            - Prefer broad cautious labels when the panel is insufficient for confident classification.
            - Include a final fallback rule only if biologically defensible.
            - Keep `features` non-exclusive; these are optional per-cell state flags, not mutually exclusive cell type assignments.
            - Keep `derived_features` dependent only on markers or feature columns defined in this YAML.
            - This is a first-pass draft for expert review, not a validated annotation schema.
            """
        ).strip()
        + "\n"
    )

    saved_path = None
    if save and output_path is None:
        output_path = (
            resolve_output_dir(adata=adata)
            / "celltype"
            / "current_prompt.txt"
        )

    metadata = adata.uns.setdefault("qxycell", {})
    stages = metadata.setdefault("stages", {})
    previous_celltyping = adata.uns.get("qxycell_celltyping", {})
    if "celltypes" not in stages and isinstance(previous_celltyping, dict):
        previous_columns = [
            previous_celltyping.get("celltype_column"),
            *previous_celltyping.get("feature_columns", []),
            *previous_celltyping.get("derived_feature_columns", []),
        ]
        stages["celltypes"] = {
            "status": "complete",
            "columns": [column for column in previous_columns if column],
            "files": [],
        }
    prepare_stage(adata, "celltype_prompt", remove_downstream_columns=False)

    if save and output_path is not None:
        output_path = Path(output_path).expanduser().resolve()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(prompt, encoding="utf-8")
        saved_path = output_path

    metadata["llm_prompt_generated"] = True
    metadata["llm_prompt_path"] = str(saved_path) if saved_path is not None else None
    complete_stage(
        adata,
        "celltype_prompt",
        files=[saved_path] if saved_path is not None else [],
        details={
            "prompt_path": str(saved_path) if saved_path is not None else None,
            "markers": markers,
        },
    )

    if print_prompt:
        if saved_path is not None:
            print(f"Saved QXYCell cell type prompt to:\n{saved_path}\n")
            print(
                "Copy the prompt below into an LLM and ask it to return YAML only. "
                "Save the returned YAML to a file such as:\n"
                f"{saved_path.parent / 'celltype_logic.yaml'}\n"
                "Then apply it with:\n"
                f"qxy.celltype(adata, {str(saved_path.parent / 'celltype_logic.yaml')!r})\n"
            )
        print(prompt)

    return prompt
