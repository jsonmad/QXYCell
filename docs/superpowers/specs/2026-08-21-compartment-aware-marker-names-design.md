# Compartment-Aware Classifier Marker Names

## Problem

QXYCell currently derives positivity-column names from classifier names alone while imported intensity variables are independently uniquified from QuPath measurement columns. For a classifier named `CD3` targeting `Nucleus: CD3 - Cy5: Median`, thresholding creates `CD3_pos`, but the source intensity is stored under a variable such as `CD3_8`. The marker heatmaps then match positivity columns to intensity variables by identical names, which selects the wrong measurement or omits the marker.

## Canonical marker naming

When QXYCell reads a classifier JSON or manual threshold definition, it will derive the public marker name from the classifier name and the compartment in its source measurement column:

| QuPath compartment | Canonical marker example |
|---|---|
| Cell or unspecified | `CD3` |
| Nucleus | `CD3-nuc` |
| Membrane | `CD3-mem` |
| Cytoplasm | `CD3-cyto` |

The statistic (`Mean` or `Median`) is not included in the public marker name. Existing classifier-name normalization remains in effect before the compartment suffix is appended.

If multiple threshold definitions still produce the same canonical marker, the existing deterministic numeric suffix convention remains responsible for uniqueness.

## Threshold data model

Threshold application will continue to preserve the original classifier name and exact `source_measurement_column`. It will additionally record, on the matched `adata.var` row, the canonical threshold marker name and corresponding positivity-column name. For example, the row whose source measurement is `Nucleus: CD3 - Cy5: Median` will link `CD3-nuc` to `CD3-nuc_pos`.

This explicit relationship is authoritative. Numeric `adata.var_names` suffixes such as `_8` remain internal uniqueness details and are not used to infer biological marker identity.

## Downstream behavior

- Thresholding creates compartment-aware positivity columns.
- `celltype_prompt()` discovers the new names from those positivity columns, so regenerated YAML examples use markers such as `CD3-nuc`.
- With `markers=None`, marker positivity and intensity heatmaps use the explicit threshold mapping and display canonical marker names.
- The intensity heatmap reads the exact intensity variable referenced by the classifier, including its configured Mean or Median statistic.
- Explicit `markers=[...]` calls retain their existing behavior.
- For older AnnData objects without the new mapping columns, plotting uses a compatibility fallback based on existing threshold metadata and classifier/source-measurement fields.

## Migration

After thresholding is rerun, cell-type YAML must reference the new compartment-aware names. A rule that previously used `CD3` for a nucleus classifier must use `CD3-nuc`. Regenerating the cell-type prompt provides the applicable names and is the preferred migration path.

Cell-compartment classifiers remain unchanged (`CD3`), limiting migration to non-cell compartments.

## Error handling

- Unknown or absent compartments receive no suffix, preserving current naming.
- Plotting ignores incomplete legacy mapping entries and uses its compatibility path.
- If no valid threshold-to-intensity mapping can be resolved, the existing clear heatmap errors remain in effect.

## Testing

Tests will verify:

1. Canonical naming for Cell, Nucleus, Membrane, Cytoplasm, and unspecified measurements.
2. Threshold application creates compartment-aware positivity columns and stores the exact mapping on the matched intensity variable.
3. Both default heatmaps include every thresholded marker, display canonical names, and use the exact configured intensity source.
4. Compatibility behavior for older AnnData objects without the new mapping fields.
5. Existing explicit marker selection remains unchanged.

## Scope

The change is limited to classifier-derived marker naming, threshold provenance, heatmap auto-selection, focused documentation, and regression tests. Imported raw intensity variable naming is not redesigned.
