# Markdown Documentation Replacement Design

## Goal

Replace the checked-in HTML documentation and GitHub Pages deployment with a
Markdown-only documentation set that renders directly in the GitHub repository.

## Documentation structure

Preserve the five current documentation entry points with Markdown equivalents:

- `docs/README.md` replaces `docs/index.html` and becomes the documentation index.
- `docs/QXYCell_overview.md` replaces `docs/QXYCell_overview.html`.
- `docs/QXYCell_flowchart.md` replaces `docs/QXYCell_flowchart.html` with a
  linear staged-workflow guide and rerun-dependency table.
- `docs/QXYCell_function_reference.md` replaces
  `docs/QXYCell_function_reference.html`.
- `docs/qxy_function_examples.md` replaces `docs/qxy_function_examples.html`.

Keep the existing focused Markdown guides. The converted pages should link to
those guides rather than duplicating large sections where a focused guide is
already authoritative.

## Content preservation

Preserve all substantive API names, parameters, command examples, workflow
stages, rerun rules, input/output descriptions, and generated example results.
Remove HTML-only navigation, CSS, layout wrappers, and decorative markup.

The generated examples HTML contains eight embedded PNG images. Extract them as
ordinary PNG files under `docs/assets/examples/` and reference them with relative
Markdown image links. Keep tables readable on GitHub; use Markdown tables where
they remain compact and fenced text blocks where very wide generated tables
would be difficult to read.

## GitHub Pages removal

Delete `.github/workflows/pages.yml` and `docs/.nojekyll`. Change the
`Documentation` project URL in `pyproject.toml` to
`https://github.com/jsonmad/QXYCell/tree/main/docs`.

Update the root README and all documentation links from `.html` targets to their
new `.md` targets. References to HTML files produced by QXYCell at runtime, such
as `dataset_summary.html`, remain unchanged because they are output-format
documentation rather than links to checked-in pages.

## Verification

- No checked-in documentation `.html` files remain under `docs/`.
- No GitHub Pages workflow, `.nojekyll`, or `jsonmad.github.io/QXYCell` project
  URL remains.
- Every repository-local Markdown link and image target resolves.
- The five Markdown replacements contain the expected major headings and API
  coverage from their HTML predecessors.
- All eight extracted PNG files are valid images and render from the examples
  document.
- `git diff --check`, Ruff, and the available package metadata checks pass.

## Scope

This is a documentation-format migration only. It does not change package
behavior, public Python APIs, runtime-generated HTML output formats, or the
existing QuPath preparation PDF.
