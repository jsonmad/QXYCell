# QXYCell Logo Spacing Design

## Objective

Use the approved balanced spacing (visual option B) between the QXYCell badge
and wordmark in every light and dark logo asset. The current PNG wordmark is too
far below the badge, while the current SVG wordmark needs slightly more space.

## Canonical geometry

The light and dark SVG files are the canonical sources. Both variants will use
the same badge and wordmark coordinates, differing only in theme colours.

At the existing 600 × 790 output size, the clear optical gap between the last
painted badge row and first painted wordmark row will target 60 px (12 SVG
viewBox units). A tolerance of 2 px is acceptable for antialiasing. Typography,
badge artwork, canvas dimensions, and horizontal alignment remain unchanged.

The wordmark placement will avoid renderer-dependent baseline behaviour. The
PNG files will be regenerated from the updated canonical SVG geometry, not
positioned or retouched independently.

## Repository-wide propagation

Update these four source assets:

- `assets/qxycell-icon.svg`
- `assets/qxycell-icon-dark.svg`
- `assets/qxycell-icon.png`
- `assets/qxycell-icon-dark.png`

Search all checked-in Markdown, HTML, PDF, generated-image, example, and source
surfaces for logo use. Rebuild every generated document that embeds one of the
changed assets, including `docs/QXYCell_QuPath_Preparation_Guide.pdf` when its
embedded logo is affected. Do not alter unrelated document content.

## Verification

- Measure the painted badge-to-wordmark gap in both PNG files at native size.
- Render both SVG files in Chrome and confirm the same balanced optical gap.
- Compare light and dark variants side by side for identical geometry.
- Confirm all four files retain their current canvas dimensions and centred
  alignment.
- Search the repository for logo references and rebuild affected generated
  documents.
- Visually inspect every page of each rebuilt PDF or other generated document.
- Run `git diff --check` and confirm only intended branding surfaces changed.

## Acceptance criteria

The four logo assets show the approved option-B spacing, light and dark variants
match geometrically, PNG and SVG presentation no longer diverges visibly, and
all checked-in logo consumers are current and visually verified.
