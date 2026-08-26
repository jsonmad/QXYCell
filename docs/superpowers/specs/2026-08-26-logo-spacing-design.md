# QXYCell Badge and Full-Lockup Design

## Objective

Use the QXYCell badge without the wordmark at the top of the GitHub README.
Retain the existing badge-plus-wordmark lockups for larger branded documents,
including guides, PDFs, and presentations.

## Asset structure

The existing light and dark full-lockup assets remain unchanged:

- `assets/qxycell-icon.svg`
- `assets/qxycell-icon-dark.svg`
- `assets/qxycell-icon.png`
- `assets/qxycell-icon-dark.png`

Add light and dark square badge variants derived from the same badge geometry:

- `assets/qxycell-badge.svg`
- `assets/qxycell-badge-dark.svg`
- `assets/qxycell-badge.png`
- `assets/qxycell-badge-dark.png`

The badge variants use a 600 x 600 canvas, preserve the existing badge artwork,
colours, border treatment, and centred alignment, and omit only the QXYCell
wordmark. The raster files are rendered from the SVG sources rather than edited
independently.

## README presentation

Replace the README's full-lockup PNG with a theme-aware HTML `picture` element.
Use the dark badge on dark colour schemes and the light badge as the fallback.
Keep the existing `# QXYCell` Markdown heading as the project name.

## Repository-wide propagation

Search all checked-in Markdown, HTML, PDF, generated-image, example, and source
surfaces for logo use. Do not replace full lockups in larger documents merely
because the README changes; the project intentionally supports both badge and
full-lockup variants. Rebuild a generated document only if its logo consumer is
changed.

## Verification

- Confirm all four new badge files are 600 x 600.
- Confirm light and dark badge SVGs contain no wordmark text.
- Render the SVG files in Chrome and compare them with their PNG counterparts.
- Render the README header in Chrome in light and dark colour schemes.
- Confirm the existing full-lockup assets are byte-for-byte unchanged.
- Search the repository for logo references and confirm each consumer uses the
  intended badge or full-lockup variant.
- Run `git diff --check` and confirm only intended branding surfaces changed.

## Acceptance criteria

The GitHub README shows a centred badge without a duplicated wordmark, switches
between light and dark badge variants with the viewer's colour scheme, and keeps
the `# QXYCell` heading. Existing full-lockup assets and larger branded documents
remain unchanged.
