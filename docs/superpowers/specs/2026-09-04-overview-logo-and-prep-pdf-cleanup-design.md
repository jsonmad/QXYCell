# Overview logo and preparation-guide cleanup

## Scope

- Render the existing QXYCell logo in `docs/QXYCell_overview.md` at a compact
  width of 96 pixels.
- Keep the existing `assets/qxycell-icon.png` source asset unchanged.
- Delete `docs/QXYCell_QuPath_Preparation_Guide.pdf`.
- Remove every published link to the deleted PDF while retaining links to
  `docs/qupath_preparation.md`.

## Implementation

Use a GitHub-compatible HTML `img` element in the Markdown overview because
standard Markdown has no portable image-size syntax. Set only `width="96"` so
the browser preserves the logo's aspect ratio.

Update the root README, documentation index, overview, function reference,
function examples, and QuPath inputs guide so preparation guidance links only
to `qupath_preparation.md`.

## Verification

- Confirm the PDF file no longer exists.
- Confirm no published document references the PDF filename.
- Confirm every local Markdown link resolves.
- Confirm only the intended documentation files changed.
