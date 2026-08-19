# QXYCell Pipeline Flowchart Redesign

## Goal

Redesign `docs/QXYCell_flowchart.html` so the workflow reads immediately as a
flowchart. The connector must be the dominant visual object, stage numbering
must match QXYCell exactly, and the Stage 3 threshold alternatives must remain
clear without complicated fork geometry.

## Root Cause

The current page uses large documentation cards connected by 2 px lines only
1.4 rem high. The nodes dominate the composition while the short connectors
blend into the card borders. This makes the page read as a vertical stack of
documentation rather than a connected process diagram. The Stage 3 fork adds
more connector geometry without improving the main route's clarity.

## Approved Visual Direction

Use a vertical pipeline route map:

- A thick teal route line runs continuously from preparation to downstream
  analysis.
- Large circular stations sit directly on the route line.
- Each station aligns with one compact information block.
- Required stations use a solid ring; optional stations use a dashed amber
  ring and a lightly tinted block; expert review uses purple.
- The route remains visually continuous through every stage, including Stage
  3.
- Text inside blocks is concise: stage/function name, one purpose sentence,
  and only the technical detail necessary to distinguish behavior.

The design remains self-contained HTML and CSS with no JavaScript or external
resources.

## Numbering

Station labels must use QXYCell's actual workflow identifiers:

1. `PREP` — prepare the QuPath project; this is not a numbered QXYCell stage.
2. A check mark — optional preflight; this is not a numbered QXYCell stage.
3. `1` — import cells.
4. `2` — add annotations.
5. `2b` — optionally remove ignored cells.
6. `3` — apply marker thresholds.
7. `4` — generate the cell-type prompt.
8. A star — expert review gate; this is not a numbered QXYCell stage.
9. `5` — assign cell types.
10. `6` — optional spatial plots.
11. An arrow — downstream analysis output.

No independent sequence such as `01`, `02`, or `03` may be introduced.

## Stage 3

Stage 3 remains one normal station and one block on the continuous main route.
The block is titled `Stage 3 · Apply marker thresholds` and says to choose
exactly one source. Inside the same block, two compact text regions describe:

- `3A · Classifier JSON` — apply QuPath classifier thresholds and save
  `thresholds/classifier_thresholds.tsv`.
- `3B · Reviewed table` — apply thresholds from one explicitly named table
  only.

The block states that 3A and 3B are alternatives and the user runs one, not
both. There is no visual fork, rejoin, branch arrow, decision diamond, or
secondary connector geometry.

## Rerun Dependencies

Keep the existing five rerun rules below the main route as a compact technical
table. The dependency section must not interrupt or cross the main connector.

## Responsive and Print Behavior

On narrow screens, the route line and stations remain visible at the left while
content blocks use the remaining width. Stage 3's two source descriptions stack
into one column. No horizontal scrolling is permitted.

For print, blocks and dependency rows avoid splitting across pages where
practical. The route may break only between blocks at a page boundary; no block
may be clipped or overlap another.

## Accessibility

Use semantic headings and ordered workflow content. State labels must remain in
text so color is never the only status indicator. Decorative route elements are
hidden from assistive technology. Preserve the skip link and visible keyboard
focus treatment for documentation links.

## Verification

Verify the redesign at desktop and narrow widths, including:

- the route line is visibly thicker than card borders;
- every station center aligns with the route line;
- stage labels match `1`, `2`, `2b`, `3`, `4`, `5`, and `6`;
- preparation, preflight, and expert review do not receive invented numbers;
- Stage 3 is one block containing both 3A and 3B descriptions;
- no legacy fork/merge connector elements remain;
- all existing local links resolve and no external resources are introduced;
- paginated print output has no split blocks, clipping, overlap, or horizontal
  overflow.

## Scope

Modify only `docs/QXYCell_flowchart.html`. Preserve the existing workflow facts,
checkpoint wording, rerun rules, and documentation links. Do not change QXYCell
runtime behavior or other documentation pages.
