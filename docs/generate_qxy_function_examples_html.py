"""Generate a self-contained HTML guide for public qxycell functions."""

from __future__ import annotations

import base64
import html
import json
import re
import shutil
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

import qxycell as qxy


ROOT = Path(__file__).resolve().parents[1]
BUILD_DIR = ROOT / "docs" / "_qxy_function_examples_build"
PROJECT_DIR = BUILD_DIR / "synthetic_qupath_export"
BASE_OUTPUT_DIR = BUILD_DIR / "outputs"
RUN_OUTPUT_DIR = BASE_OUTPUT_DIR / "synthetic_run_000000_0000"
WORKFLOW_OUTPUT_DIR = BASE_OUTPUT_DIR / "synthetic_workflow_run_000000_0001"
HTML_PATH = ROOT / "docs" / "qxy_function_examples.html"


def reset_build() -> None:
    if BUILD_DIR.exists():
        shutil.rmtree(BUILD_DIR)
    PROJECT_DIR.mkdir(parents=True)
    BASE_OUTPUT_DIR.mkdir(parents=True)


def classifier_json(measurement: str, threshold: float) -> dict:
    return {
        "function": {
            "classifier_fun": "ClassifyByMeasurementFunction",
            "measurement": measurement,
            "threshold": threshold,
        }
    }


def polygon_feature(object_type: str, name: str, coords: list[list[float]]) -> dict:
    return {
        "type": "Feature",
        "properties": {"objectType": object_type, "name": name},
        "geometry": {"type": "Polygon", "coordinates": [coords]},
    }


def cell_feature(object_id: str, x: float, y: float, half_width: float = 1.5) -> dict:
    coords = [
        [x - half_width, y - half_width],
        [x + half_width, y - half_width],
        [x + half_width, y + half_width],
        [x - half_width, y + half_width],
        [x - half_width, y - half_width],
    ]
    feature = polygon_feature("cell", "", coords)
    feature["id"] = object_id
    return feature


def write_geojson(path: Path, features: list[dict]) -> None:
    path.write_text(
        json.dumps({"type": "FeatureCollection", "features": features}, indent=2),
        encoding="utf-8",
    )


def make_synthetic_project() -> None:
    rng = np.random.default_rng(7)
    rows = []
    cell_features_by_image: dict[str, list[dict]] = {"img_a": [], "img_b": []}
    images = [
        ("img_a.ome.tiff", "Sample-A", 0.0),
        ("img_b.ome.tiff", "Sample-B", 120.0),
    ]
    for image_name, sample_name, x_offset in images:
        image_key = Path(image_name).stem.removesuffix(".ome")
        for i in range(36):
            object_id = f"{image_key}_cell_{i:03d}"
            cluster = i % 4
            cx = x_offset + 10 + (i % 6) * 10 + rng.normal(0, 1.0)
            cy = 10 + (i // 6) * 10 + rng.normal(0, 1.0)
            cd3 = 1.0 if cluster in {0, 1} else 0.0
            cd8 = 1.0 if cluster == 0 else 0.0
            panck = 1.0 if cluster == 2 else 0.0
            cd68 = 1.0 if cluster == 3 else 0.0
            rows.append(
                {
                    "Image": image_name,
                    "Object ID": object_id,
                    "Centroid X um": cx,
                    "Centroid Y um": cy,
                    "Centroid X \u00b5m": cx,
                    "Centroid Y \u00b5m": cy,
                    "TMA Core": f"{sample_name}-Core-{(i % 2) + 1}",
                    "CD3: Mean": cd3 + rng.normal(0.1, 0.04),
                    "CD8: Mean": cd8 + rng.normal(0.1, 0.04),
                    "PanCK: Mean": panck + rng.normal(0.1, 0.04),
                    "CD68: Mean": cd68 + rng.normal(0.1, 0.04),
                }
            )
            cell_features_by_image[image_key].append(cell_feature(object_id, cx, cy))

    pd.DataFrame(rows).to_csv(PROJECT_DIR / "detections.tsv", sep="\t", index=False)
    for marker in ("CD3", "CD8", "PanCK", "CD68"):
        (PROJECT_DIR / f"{marker}.json").write_text(
            json.dumps(classifier_json(f"{marker}: Mean", 0.5), indent=2),
            encoding="utf-8",
        )

    write_geojson(
        PROJECT_DIR / "img_a.geojson",
        [
            polygon_feature("annotation", "Sample-A", [[0, 0], [75, 0], [75, 75], [0, 75], [0, 0]]),
            polygon_feature("annotation", "Immune zone", [[0, 0], [45, 0], [45, 75], [0, 75], [0, 0]]),
            polygon_feature("annotation", "Ignore", [[5, 5], [18, 5], [18, 18], [5, 18], [5, 5]]),
            polygon_feature("tmaCore", "A_core_1", [[0, 0], [38, 0], [38, 75], [0, 75], [0, 0]]),
            polygon_feature("tmaCore", "A_core_2", [[38, 0], [75, 0], [75, 75], [38, 75], [38, 0]]),
        ],
    )
    write_geojson(
        PROJECT_DIR / "img_b.geojson",
        [
            polygon_feature("annotation", "Sample-B", [[115, 0], [190, 0], [190, 75], [115, 75], [115, 0]]),
            polygon_feature("annotation", "Tumour bed", [[145, 0], [190, 0], [190, 75], [145, 75], [145, 0]]),
            polygon_feature("tmaCore", "B_core_1", [[115, 0], [153, 0], [153, 75], [115, 75], [115, 0]]),
            polygon_feature("tmaCore", "B_core_2", [[153, 0], [190, 0], [190, 75], [153, 75], [153, 0]]),
        ],
    )
    for image_key, features in cell_features_by_image.items():
        write_geojson(PROJECT_DIR / f"{image_key}-cells.geojson", features)


def data_uri(path: Path) -> str:
    suffix = path.suffix.lower()
    mime = "image/png" if suffix == ".png" else "image/svg+xml"
    payload = base64.b64encode(path.read_bytes()).decode("ascii")
    return f"data:{mime};base64,{payload}"


def code_block(code: str) -> str:
    return f"<pre><code>{html.escape(textwrap.dedent(code).strip())}</code></pre>"


def table_html(df: pd.DataFrame, max_rows: int = 8) -> str:
    return df.head(max_rows).to_html(index=False, border=0, classes="data-table")


def compact_dict(d: dict, keys: list[str]) -> pd.DataFrame:
    return pd.DataFrame([{"field": key, "value": d.get(key)} for key in keys])


def display_path(path: str | Path) -> str:
    """Return a repository-relative path suitable for public generated docs."""

    resolved = Path(path).resolve()
    try:
        return f"REPOSITORY_ROOT/{resolved.relative_to(ROOT).as_posix()}"
    except ValueError:
        return resolved.name


def build_examples() -> dict:
    threshold_path = qxy.generate_threshold_table(
        PROJECT_DIR,
        output_dir=RUN_OUTPUT_DIR / "threshold_example",
    )
    report = qxy.check(PROJECT_DIR, output_dir=RUN_OUTPUT_DIR, count_rows=True)
    logic = {
        "rules": [
            {"name": "CD8 T", "positive": ["CD3", "CD8"]},
            {"name": "Other T", "positive": ["CD3"]},
            {"name": "Tumour", "positive": ["PanCK"]},
            {"name": "Macrophage", "positive": ["CD68"]},
        ],
        "features": {"immune_marker_positive": {"any_of": ["CD3", "CD68"]}},
        "derived_features": {},
    }
    logic_path = BUILD_DIR / "celltype_logic.yaml"
    logic_path.write_text(yaml.safe_dump(logic, sort_keys=False), encoding="utf-8")

    adata = qxy.run(
        PROJECT_DIR,
        output_dir=RUN_OUTPUT_DIR,
        pixel_size_um=1.0,
        verbose=False,
    )
    threshold_summary = qxy.threshold(
        adata,
        PROJECT_DIR,
        output_dir=RUN_OUTPUT_DIR,
        verbose=False,
    )
    apply_thresholds_summary = qxy.apply_thresholds(
        adata.copy(),
        PROJECT_DIR,
        output_dir=RUN_OUTPUT_DIR,
        verbose=False,
    )
    celltype_logic_copy = RUN_OUTPUT_DIR / "celltype" / "celltype_logic.yaml"
    celltype_logic_copy.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(logic_path, celltype_logic_copy)
    loaded_logic = qxy.load_celltype_logic(logic_path)
    latest_logic_path = qxy.find_latest_celltype_yaml(RUN_OUTPUT_DIR / "celltype")
    prompt_path = BUILD_DIR / "celltype_prompt.txt"
    prompt = qxy.celltype_prompt(
        adata,
        context="Small synthetic QXYCell demo with CD3, CD8, PanCK, and CD68.",
        output_path=prompt_path,
        print_prompt=False,
    )

    metadata = pd.DataFrame(
        {
            "Sample": ["Sample-A", "Sample-B"],
            "condition": ["treated", "control"],
            "batch": ["B1", "B1"],
        }
    )
    metadata_summary = qxy.add_metadata(
        adata,
        metadata,
        sample_col="Sample",
        output_dir=RUN_OUTPUT_DIR,
        verbose=False,
    )
    clean = qxy.remove_ignore(adata, copy=True, verbose=False)
    artifact_clean = qxy.remove_annotations(
        adata,
        text="Immune",
        copy=True,
        verbose=False,
    )
    sample_adata = adata.copy()
    sample_adata.obs["annotation__Sample_A"] = sample_adata.obs["Sample"].astype(str).eq("Sample-A")
    sample_adata.obs["annotation__Sample_B"] = sample_adata.obs["Sample"].astype(str).eq("Sample-B")
    sample_adata.uns["qxycell_annotation_labels"] = {
        **sample_adata.uns.get("qxycell_annotation_labels", {}),
        "annotation__Sample_A": "Sample-A",
        "annotation__Sample_B": "Sample-B",
    }
    sample_summary = qxy.assign_samples(sample_adata, overwrite=True, verbose=False)
    annotation_adata = adata.copy()
    annotation_summary = qxy.assign_annotations(
        annotation_adata,
        ["Immune zone", "Tumour bed"],
        target_col="Region",
        drop=False,
        verbose=False,
    )
    coreid_summary = qxy.assign_core_ids_from_measurements(adata, verbose=False)
    n_polygons = qxy.load_cell_polygons(
        adata,
        PROJECT_DIR,
        pixel_size_um=1.0,
        verbose=False,
    )

    celltype_summary = qxy.celltype(
        adata,
        loaded_logic,
        celltype_dir=RUN_OUTPUT_DIR / "celltype_rerun",
        verbose=False,
    )
    apply_celltypes_summary = qxy.apply_celltypes(
        adata,
        loaded_logic,
        celltype_dir=RUN_OUTPUT_DIR / "celltype_apply_alias",
        verbose=False,
    )
    qxy.cn_knn(adata, k=5)
    qxy.cn_kmeans(adata, n_cn=3, random_state=1)
    cn_labels = qxy.cn_name(adata, output_dir=RUN_OUTPUT_DIR / "cn", verbose=False)
    qc_summary = qxy.qc(adata, sample_col="Sample", output_dir=RUN_OUTPUT_DIR, verbose=False)

    plot_paths = {}
    plot_paths["stacked_bar"] = qxy.plot_stacked_bar(
        adata,
        sample_col="Sample",
        output_dir=RUN_OUTPUT_DIR,
        show_axis_labels=False,
        show=False,
        verbose=False,
        save_png=True,
        save_pdf=False,
    )
    plot_paths["spatial"] = qxy.plot_spatial(
        adata,
        sample_col="Image",
        output_dir=RUN_OUTPUT_DIR,
        show=False,
        verbose=False,
        point_size=16,
        scale_bar=False,
        auto_figsize=True,
        save_pdf=False,
    )
    plot_paths["cell_boundaries"] = qxy.plot_cell_boundaries(
        adata,
        sample_col="Sample",
        output_dir=RUN_OUTPUT_DIR,
        show=False,
        verbose=False,
        scale_bar=False,
        auto_figsize=True,
        save_pdf=False,
        label_celltypes="Tumour",
        label_color="black",
        label_offset_um=(6.0, 6.0),
    )
    plot_paths["annotation_polygons"] = qxy.plot_annotation_polygons(
        adata,
        project_dir=PROJECT_DIR,
        output_dir=RUN_OUTPUT_DIR,
        pixel_size_um=1.0,
        show=False,
        verbose=False,
    )
    plot_paths["marker_positivity_heatmap"] = qxy.plot_marker_positivity_heatmap(
        adata,
        output_dir=RUN_OUTPUT_DIR,
        show=False,
        verbose=False,
        annotate=True,
        save_png=True,
        save_pdf=False,
    )
    plot_paths["marker_intensity_heatmap"] = qxy.plot_marker_intensity_heatmap(
        adata,
        output_dir=RUN_OUTPUT_DIR,
        show=False,
        verbose=False,
        annotate=True,
        save_png=True,
        save_pdf=False,
    )
    plot_paths["marker_heatmap_compatibility"] = qxy.plot_marker_heatmap(
        adata,
        values="positivity",
        output_dir=RUN_OUTPUT_DIR,
        show=False,
        verbose=False,
        save_png=True,
        save_pdf=False,
    )
    plot_paths["cn_heatmap"] = qxy.plot_cn_heatmap(
        adata,
        sample_col="Sample",
        condition_col="condition",
        output_dir=RUN_OUTPUT_DIR,
        show=False,
        verbose=False,
        annotate=True,
        save_png=True,
        save_pdf=False,
    )

    saved_path = qxy.save(adata, output_dir=RUN_OUTPUT_DIR, verbose=False)
    loaded = qxy.load(RUN_OUTPUT_DIR)
    latest = qxy.load_latest(BASE_OUTPUT_DIR)
    workflow_adata = qxy.workflow(
        PROJECT_DIR,
        output_dir=WORKFLOW_OUTPUT_DIR,
        sample_metadata=metadata,
        sample_col="Sample",
        celltype_logic=logic_path,
        remove_ignore_cells=True,
        make_qc=True,
        make_plots=True,
        plot_sample_col="Sample",
        verbose=False,
    )

    return {
        "report": report,
        "threshold_path": threshold_path,
        "threshold_summary": threshold_summary,
        "apply_thresholds_summary": apply_thresholds_summary,
        "adata": adata,
        "clean": clean,
        "artifact_clean": artifact_clean,
        "workflow_adata": workflow_adata,
        "logic_path": logic_path,
        "latest_logic_path": latest_logic_path,
        "prompt": prompt,
        "metadata_summary": metadata_summary,
        "sample_summary": sample_summary,
        "annotation_summary": annotation_summary,
        "coreid_summary": coreid_summary,
        "n_polygons": n_polygons,
        "celltype_summary": celltype_summary,
        "apply_celltypes_summary": apply_celltypes_summary,
        "cn_labels": cn_labels,
        "qc_summary": qc_summary,
        "plot_paths": plot_paths,
        "saved_path": saved_path,
        "loaded_shape": loaded.shape,
        "latest_shape": latest.shape,
    }


def image_from_glob(pattern: str) -> str:
    matches = sorted(RUN_OUTPUT_DIR.rglob(pattern))
    if not matches:
        return ""
    return f'<img class="figure" src="{data_uri(matches[0])}" alt="{html.escape(matches[0].name)}">'


def make_html(examples: dict) -> str:
    adata = examples["adata"]
    report = examples["report"]
    qc_summary = examples["qc_summary"]
    qc_overview = dict(
        zip(
            qc_summary["tables"]["overview"]["metric"],
            qc_summary["tables"]["overview"]["value"],
            strict=True,
        )
    )
    obs_cols = ["Image", "Object ID", "Sample", "CoreID", "celltype", "cn", "condition"]
    obs_preview = adata.obs[[c for c in obs_cols if c in adata.obs.columns]].head(8).reset_index()
    marker_preview = adata.obs[
        [c for c in adata.obs.columns if c.endswith("_pos")]
    ].head(8).reset_index()
    functions = [
        ("qxy.run()", "Import QuPath measurements, spatial coordinates, annotation labels, cell polygons, and an optional TMA Core column. Thresholding and cell typing are off by default; the run log records whether each was applied and its source.", "adata = qxy.run(project_dir, output_dir=out_dir, pixel_size_um=1.0)", table_html(obs_preview)),
        ("qxy.threshold()", "Apply threshold definitions to an imported AnnData object, add <marker>_pos columns, and rename previous cell type labels with a __stale_celltype suffix.", "summary = qxy.threshold(adata, project_dir, output_dir=out_dir)", table_html(compact_dict(examples["threshold_summary"], ["threshold_source", "n_threshold_definitions", "n_pos_columns"]))),
        ("qxy.apply_thresholds()", "Long-name alias for qxy.threshold().", "summary = qxy.apply_thresholds(adata, project_dir, output_dir=out_dir)", table_html(compact_dict(examples["apply_thresholds_summary"], ["threshold_source", "n_threshold_definitions", "n_pos_columns"]))),
        ("qxy.workflow()", "Run the common Python workflow in one call.", "adata = qxy.workflow(project_dir, sample_metadata=metadata, celltype_logic='celltype_logic.yaml')", f"<p>Workflow result: {examples['workflow_adata'].n_obs} cells after Ignore removal.</p>"),
        (
            "qxy.check()",
            "Inspect the export folder, list every annotation and planned AnnData assignment, and audit competing classifier thresholds without silently choosing one. Check never applies thresholds or cell typing and never generates an LLM prompt.",
            "report = qxy.check(project_dir, output_dir=out_dir, count_rows=True)",
            table_html(
                pd.DataFrame(
                    [
                        {"metric": "ok", "value": report.ok},
                        {
                            "metric": "threshold source",
                            "value": report.active_threshold_source
                            or report.active_threshold_source_kind,
                        },
                        {"metric": "measurement files", "value": len(report.measurement_files)},
                        {"metric": "geojson files", "value": len(report.geojson_files)},
                    ]
                )
            ),
        ),
        (
            "qxy.generate_threshold_table()",
            "Create a fresh timestamped threshold table. Conflicted channels retain candidate provenance and blank per-image cells until reviewed.",
            "threshold_path = qxy.generate_threshold_table(project_dir, output_dir=out_dir)",
            f"<p><code>{html.escape(display_path(examples['threshold_path']))}</code></p>",
        ),
        ("qxy.qc()", "Write QC tables for markers, annotations, samples, and cell types.", "qc_summary = qxy.qc(adata, sample_col='Sample', output_dir=out_dir)", table_html(compact_dict({"output_dir": display_path(qc_summary["output_dir"]), **qc_overview}, ["output_dir", "n_cells", "sample_col"]))),
        ("qxy.save()", "Save an AnnData object to H5AD.", "h5ad_path = qxy.save(adata, output_dir=out_dir)", f"<p><code>{html.escape(display_path(examples['saved_path']))}</code></p>"),
        ("qxy.load()", "Load an H5AD from a file or QXYCell output folder.", "adata = qxy.load(out_dir)", f"<p>Loaded shape: <code>{examples['loaded_shape']}</code></p>"),
        ("qxy.load_latest()", "Load the newest qxy_outputs_* or <project>_run_* folder in a base directory.", "adata = qxy.load_latest(base_output_dir)", f"<p>Latest shape: <code>{examples['latest_shape']}</code></p>"),
        ("qxy.assign_samples()", "Create one Sample column from sample annotation labels.", "summary = qxy.assign_samples(adata)", table_html(compact_dict(examples["sample_summary"], ["sample_col", "n_assigned_cells", "n_conflicting_cells"]))),
        ("qxy.assign_annotations()", "Collapse selected boolean annotation columns into one categorical observation column.", "summary = qxy.assign_annotations(adata, ['Immune zone', 'Tumour bed'], target_col='Region', drop=False)", table_html(compact_dict(examples["annotation_summary"], ["target_col", "n_assigned_cells", "n_unassigned_cells", "n_conflicting_cells"]))),
        ("qxy.assign_core_ids_from_measurements()", "Create CoreID only from QuPath's measurement-table TMA Core column. Parent and annotations are never used.", "summary = qxy.assign_core_ids_from_measurements(adata)", table_html(compact_dict(examples["coreid_summary"], ["target_col", "available_source_cols", "n_assigned_cells", "n_unassigned_cells"]))),
        ("qxy.remove_ignore()", "Remove cells in annotation columns containing Ignore.", "clean = qxy.remove_ignore(adata, copy=True)", f"<p>Before: {adata.n_obs} cells. After copy: {examples['clean'].n_obs} cells.</p>"),
        ("qxy.remove_annotations()", "Remove cells in annotation columns matching a custom text string.", "clean = qxy.remove_annotations(adata, text='Immune', copy=True)", f"<p>Before: {adata.n_obs} cells. After copy: {examples['artifact_clean'].n_obs} cells.</p>"),
        ("qxy.load_cell_polygons()", "Load cell boundary polygons from GeoJSON into cell_polygon_wkt.", "n = qxy.load_cell_polygons(adata, project_dir, pixel_size_um=1.0)", f"<p>Matched polygons: <code>{examples['n_polygons']}</code></p>"),
        ("qxy.add_metadata()", "Attach sample-level metadata to cells.", "qxy.add_metadata(adata, metadata, sample_col='Sample')", table_html(compact_dict(examples["metadata_summary"], ["n_matched_samples", "added_columns"]))),
        ("qxy.celltype()", "Apply ordered cell type rules.", "summary = qxy.celltype(adata, 'celltype_logic.yaml')", table_html(compact_dict(examples["celltype_summary"], ["n_rules", "unknown_count", "celltype_column"]))),
        ("qxy.apply_celltypes()", "Public long-name alias for applying ordered cell type rules.", "summary = qxy.apply_celltypes(adata, 'celltype_logic.yaml')", table_html(compact_dict(examples["apply_celltypes_summary"], ["n_rules", "unknown_count", "celltype_column"]))),
        ("qxy.load_celltype_logic()", "Load cell typing YAML as a Python dict.", "logic = qxy.load_celltype_logic('celltype_logic.yaml')", code_block(yaml.safe_dump(qxy.load_celltype_logic(examples["logic_path"]), sort_keys=False))),
        ("qxy.find_latest_celltype_yaml()", "Find the newest saved cell type YAML.", "path = qxy.find_latest_celltype_yaml(out_dir / 'celltype')", f"<p><code>{html.escape(display_path(examples['latest_logic_path']))}</code></p>"),
        ("qxy.CheckReport", "Result type returned by qxy.check().", "isinstance(report, qxy.CheckReport)", f"<p><code>{isinstance(report, qxy.CheckReport)}</code>; errors: <code>{report.n_errors}</code>; warnings: <code>{report.n_warnings}</code></p>"),
        ("qxy.celltype_prompt()", "Generate a prompt for drafting cell type rules.", "prompt = qxy.celltype_prompt(adata, print_prompt=False)", code_block(examples["prompt"][:900] + "\\n...")),
        ("qxy.cn_knn()", "Build cellular neighbourhood composition profiles.", "qxy.cn_knn(adata, k=5)", f"<p><code>adata.obsm['cn_profile']</code> shape: {adata.obsm['cn_profile'].shape}</p>"),
        ("qxy.cn_kmeans()", "Cluster neighbourhood profiles into CNs.", "qxy.cn_kmeans(adata, n_cn=3)", table_html(adata.obs["cn"].value_counts().rename_axis("cn").reset_index(name="n_cells"))),
        ("qxy.cn_name()", "Rename CN clusters by their cell type composition.", "labels = qxy.cn_name(adata)", table_html(examples["cn_labels"])),
        ("qxy.plot_stacked_bar()", "Plot cell type frequencies by sample.", "qxy.plot_stacked_bar(adata, sample_col='Sample', show_axis_labels=False, show=False)", image_from_glob("*stacked_bar*.png")),
        ("qxy.plot_spatial()", "Plot spatial cell locations by category. Missing sample labels are excluded by default; underlay_adata can provide a full-data grey underlay for a filtered plot object.", "qxy.plot_spatial(adata_cn, underlay_adata=adata, sample_col='Sample', category_col='cn', show=False)", image_from_glob("*spatial*.png")),
        ("qxy.plot_cell_boundaries()", "Plot cell boundary polygons by category.", "qxy.plot_cell_boundaries(adata, sample_col='Sample', label_celltypes='Tumour', save_pdf=False, show=False)", image_from_glob("*cell_boundaries*.png")),
        ("qxy.plot_annotation_polygons()", "Reload QuPath annotation geometry for PNG-only QC. Defaults to a low-resolution cell-density underlay and boundary-only polygons (fill=False).", "qxy.plot_annotation_polygons(adata, show=False)", image_from_glob("*annotation_polygons*.png")),
        ("qxy.plot_marker_positivity_heatmap()", "Plot positive-cell count divided by category cell count for each thresholded marker (0–1).", "qxy.plot_marker_positivity_heatmap(adata, annotate=True, show=False)", image_from_glob("*marker_heatmap_positivity.png")),
        ("qxy.plot_marker_intensity_heatmap()", "Plot category mean intensity followed by per-marker Z-scoring. This is not a median or count heatmap.", "qxy.plot_marker_intensity_heatmap(adata, annotate=True, show=False)", image_from_glob("*marker_heatmap_intensity.png")),
        ("qxy.plot_marker_heatmap()", "Compatibility entry point for older code. New code should call the explicit positivity or intensity function.", "qxy.plot_marker_heatmap(adata, values='positivity', show=False)", image_from_glob("*marker_heatmap_positivity.png")),
        ("qxy.plot_cn_heatmap()", "Plot CN abundance heatmap across samples.", "qxy.plot_cn_heatmap(adata, sample_col='Sample', condition_col='condition', show=False)", image_from_glob("*cn_heatmap*.png")),
    ]
    cards = []
    for name, purpose, code, output in functions:
        cards.append(
            f"""
            <section class="card">
              <h3>{html.escape(name)}</h3>
              <p>{html.escape(purpose)}</p>
              {code_block(code)}
              <div class="output">{output}</div>
            </section>
            """.strip()
        )

    generated_html = f"""<!doctype html>
<!-- Generated by docs/generate_qxy_function_examples_html.py from synthetic data. Do not edit directly. -->
<html lang="en">
<head>
<meta charset="utf-8">
<title>QXYCell Function Examples</title>
<style>
  :root {{
    --ink: #18212b;
    --muted: #5d6876;
    --line: #d8dee7;
    --panel: #f7f9fb;
    --accent: #1b6b6f;
    --accent-soft: #e3f0ef;
  }}
  body {{
    margin: 0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Arial, sans-serif;
    color: var(--ink);
    background: white;
    line-height: 1.45;
  }}
  main {{
    max-width: 1100px;
    margin: 0 auto;
    padding: 40px 28px 72px;
  }}
  header {{
    border-bottom: 3px solid var(--accent);
    padding-bottom: 18px;
    margin-bottom: 26px;
  }}
  h1 {{
    margin: 0 0 8px;
    font-size: 34px;
    letter-spacing: 0;
  }}
  h2 {{
    margin-top: 34px;
    border-bottom: 1px solid var(--line);
    padding-bottom: 6px;
    font-size: 22px;
  }}
  h3 {{
    margin: 0 0 6px;
    color: var(--accent);
    font-size: 18px;
  }}
  p {{
    margin: 0 0 10px;
  }}
  .lead {{
    color: var(--muted);
    font-size: 16px;
    max-width: 900px;
  }}
  .summary {{
    display: grid;
    grid-template-columns: repeat(4, minmax(0, 1fr));
    gap: 10px;
    margin: 22px 0;
  }}
  .metric {{
    background: var(--accent-soft);
    border: 1px solid #c6dfdd;
    padding: 12px;
    border-radius: 6px;
  }}
  .metric strong {{
    display: block;
    font-size: 20px;
  }}
  .grid {{
    display: grid;
    grid-template-columns: repeat(2, minmax(0, 1fr));
    gap: 14px;
  }}
  .card {{
    border: 1px solid var(--line);
    border-radius: 8px;
    padding: 16px;
    background: #fff;
    break-inside: avoid;
  }}
  pre {{
    background: #111923;
    color: #e8eef5;
    padding: 12px;
    border-radius: 6px;
    overflow-x: auto;
    font-size: 12px;
    line-height: 1.4;
  }}
  code {{
    font-family: "SFMono-Regular", Consolas, monospace;
  }}
  .output {{
    margin-top: 10px;
    background: var(--panel);
    border: 1px solid var(--line);
    border-radius: 6px;
    padding: 10px;
    overflow-x: auto;
  }}
  .data-table {{
    border-collapse: collapse;
    width: 100%;
    font-size: 12px;
  }}
  .data-table th {{
    text-align: left;
    background: #edf2f7;
    border-bottom: 1px solid var(--line);
    padding: 6px;
  }}
  .data-table td {{
    border-bottom: 1px solid #e7ebf0;
    padding: 6px;
    vertical-align: top;
  }}
  .figure {{
    display: block;
    max-width: 100%;
    height: auto;
    background: white;
  }}
  .note {{
    border-left: 4px solid var(--accent);
    padding: 10px 12px;
    background: var(--accent-soft);
    margin: 14px 0;
  }}
  @media print {{
    main {{ padding: 20px; }}
    .grid {{ grid-template-columns: 1fr; }}
    .card {{ page-break-inside: avoid; }}
  }}
</style>
</head>
<body>
<main>
<header>
  <h1>QXYCell Function Examples</h1>
  <p class="lead">A single-file HTML reference for the public <code>qxy</code> API. Examples were generated from a small synthetic QuPath export and include real tables and plots produced by the current package code.</p>
</header>
<section class="summary">
  <div class="metric"><span>Cells</span><strong>{adata.n_obs}</strong></div>
  <div class="metric"><span>Markers</span><strong>{adata.n_vars}</strong></div>
  <div class="metric"><span>Samples</span><strong>{adata.obs['Sample'].nunique()}</strong></div>
  <div class="metric"><span>Functions</span><strong>{len(functions)}</strong></div>
</section>
<section class="note">
  <p><strong>Annotation rule shown here:</strong> sample annotations collapse into one <code>Sample</code> column; Ignore and other annotations remain boolean <code>annotation__*</code> columns; only the measurement-table <code>TMA Core</code> column creates <code>CoreID</code>. GeoJSON TMA core objects are reported but are not used for core assignment.</p>
</section>
<h2>Example Dataset Preview</h2>
{table_html(obs_preview)}
<h2>Marker Positivity Preview</h2>
{table_html(marker_preview)}
<h2>Function Examples</h2>
<div class="grid">
{''.join(cards)}
</div>
</main>
</body>
</html>
"""
    public_html = generated_html.replace(str(ROOT), "REPOSITORY_ROOT")
    return re.sub(r"thresholds_\d{6}-\d{4}", "thresholds_YYMMDD-HHMM", public_html)


def main() -> None:
    reset_build()
    make_synthetic_project()
    examples = build_examples()
    HTML_PATH.write_text(make_html(examples), encoding="utf-8")
    print(HTML_PATH)


if __name__ == "__main__":
    main()
