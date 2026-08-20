# QXYCell flowchart tool-color design

## Goal

Make the workflow's execution context immediately visible without adding visual clutter. PREP must read as a QuPath activity, the main numbered stages must identify their tools, and status labels must be limited to meaningful exceptions.

## Visual system

- Give PREP a distinct blue QuPath accent across its station, card border, and `QuPath` badge.
- Keep QXYCell Python stages teal and label them `Python · QXYCell`.
- Give Stage 4 an indigo `Python · QXYCell + LLM` badge.
- Keep the domain-review step purple and label it `Domain expert · YAML`.
- Keep optional stages amber for status while also showing their teal `Python · QXYCell` tool badge.
- Keep Stage 3's `Choose one route` badge alongside its Python tool badge.
- Keep the downstream status badge and add a `Python analysis` tool badge.

## Label changes

- Remove every `Required` badge from workflow cards.
- Remove the page-level `Required workflow` badge and the `Required` legend entry.
- Retain `Optional`, `Choose one route`, `Domain expert review`, and `Downstream` where they communicate status or responsibility.

## Implementation boundaries

Change only `docs/QXYCell_flowchart.html`. Reuse the existing badge, station, and card structures with small CSS additions rather than introducing a new component system. Do not alter the agreed workflow wording, step order, connectors, rerun dependencies, or print structure.

## Verification

- Confirm no workflow badge contains `Required`.
- Confirm every workflow card has the agreed tool/context badge.
- Confirm PREP has its own blue station and card accent.
- Check desktop and mobile layouts for overflow, connector alignment, and readable badge wrapping.
- Confirm the A4 print layout remains usable.
