# Markdown Documentation Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace all checked-in HTML documentation and GitHub Pages deployment with equivalent GitHub-rendered Markdown documentation.

**Architecture:** Convert each HTML entry point to a same-purpose Markdown document, with `docs/README.md` as the index. Extract the eight embedded example figures into normal PNG assets, then remove Pages infrastructure and update all repository navigation to Markdown targets.

**Tech Stack:** Markdown, Pandoc, Python standard library, Pillow, Git

## Global Constraints

- Preserve substantive API names, parameters, commands, workflow stages, rerun rules, input/output descriptions, and generated example results.
- Do not change package behavior, public Python APIs, runtime-generated HTML output formats, or the QuPath preparation PDF.
- Keep existing focused Markdown guides authoritative and avoid new documentation tooling or dependencies.
- Preserve user-authored repository changes and keep edits limited to documentation and metadata links.

---

### Task 1: Convert the four textual HTML pages

**Files:**
- Create: `docs/README.md`
- Create: `docs/QXYCell_overview.md`
- Create: `docs/QXYCell_flowchart.md`
- Create: `docs/QXYCell_function_reference.md`
- Delete: `docs/index.html`
- Delete: `docs/QXYCell_overview.html`
- Delete: `docs/QXYCell_flowchart.html`
- Delete: `docs/QXYCell_function_reference.html`

**Interfaces:**
- Consumes: the headings, prose, tables, lists, and code blocks in the four HTML source pages.
- Produces: four relative-link Markdown documents rendered directly by GitHub.

- [ ] **Step 1: Record source coverage**

Run:

```bash
rg -o '<h[1-6][^>]*>.*</h[1-6]>' docs/index.html docs/QXYCell_overview.html docs/QXYCell_flowchart.html docs/QXYCell_function_reference.html
```

Expected: headings for the index, overview, staged workflow, function groups, commands, and outputs.

- [ ] **Step 2: Convert the documents with Pandoc**

Run one conversion per file with wrapping disabled:

```bash
pandoc --from=html --to=gfm --wrap=none docs/index.html -o docs/README.md
pandoc --from=html --to=gfm --wrap=none docs/QXYCell_overview.html -o docs/QXYCell_overview.md
pandoc --from=html --to=gfm --wrap=none docs/QXYCell_flowchart.html -o docs/QXYCell_flowchart.md
pandoc --from=html --to=gfm --wrap=none docs/QXYCell_function_reference.html -o docs/QXYCell_function_reference.md
```

Expected: four Markdown files containing the source headings and content without CSS.

- [ ] **Step 3: Normalize Markdown structure and relative links**

Use direct edits to remove Pandoc wrapper artifacts, normalize headings and tables, retain language-tagged code fences, and replace internal `.html` targets with `.md` targets.

- [ ] **Step 4: Delete the replaced HTML files**

Delete the four source HTML files only after their content coverage has been checked against the Markdown outputs.

- [ ] **Step 5: Verify task output**

Run:

```bash
rg -n '^#{1,6} ' docs/README.md docs/QXYCell_overview.md docs/QXYCell_flowchart.md docs/QXYCell_function_reference.md
rg -n 'qxy\.(check|import_cells|add_annotations|threshold_from_classifiers|threshold_from_table|celltype_prompt|celltype|plot_spatial)' docs/QXYCell_overview.md docs/QXYCell_flowchart.md docs/QXYCell_function_reference.md
```

Expected: all major sections and core API names remain present.

### Task 2: Convert generated function examples and extract figures

**Files:**
- Create: `docs/qxy_function_examples.md`
- Create: `docs/assets/examples/*.png`
- Delete: `docs/qxy_function_examples.html`

**Interfaces:**
- Consumes: generated prose, code, tables, output blocks, and eight `data:image/png;base64,...` figures.
- Produces: a readable Markdown examples reference and eight relative PNG assets.

- [ ] **Step 1: Extract embedded figures**

Use a temporary Python script that parses each embedded PNG data URI, decodes it, derives a safe filename from the image `alt` value, and writes it to `docs/assets/examples/`. Reject duplicate filenames and non-PNG payloads.

- [ ] **Step 2: Replace data URIs with relative image paths**

Create a temporary HTML copy in `/private/tmp`, replace each embedded `src` with `assets/examples/<filename>`, and retain its `alt` description.

- [ ] **Step 3: Convert the temporary HTML to Markdown**

Run:

```bash
pandoc --from=html --to=gfm --wrap=none /private/tmp/qxy_function_examples_externalized.html -o docs/qxy_function_examples.md
```

Expected: generated example headings, calls, outputs, tables, and eight Markdown image references.

- [ ] **Step 4: Normalize large generated tables**

Keep compact tables as Markdown tables. Convert very wide or multiline tables to fenced `text` blocks so GitHub renders them predictably without changing values.

- [ ] **Step 5: Delete the generated HTML source**

Delete `docs/qxy_function_examples.html` after verifying the Markdown and image assets.

- [ ] **Step 6: Verify task output**

Run:

```bash
test "$(find docs/assets/examples -type f -name '*.png' | wc -l | tr -d ' ')" = 8
rg -n '^### qxy\.' docs/qxy_function_examples.md
rg -n '!\[[^]]*\]\(assets/examples/[^)]+\.png\)' docs/qxy_function_examples.md
```

Expected: eight PNG files, complete function headings, and eight relative image references.

### Task 3: Remove GitHub Pages and update repository navigation

**Files:**
- Delete: `.github/workflows/pages.yml`
- Delete: `docs/.nojekyll`
- Modify: `README.md`
- Modify: `docs/plotting.md`
- Modify: `docs/qupath_inputs.md`
- Modify: `pyproject.toml`

**Interfaces:**
- Consumes: the Markdown filenames created in Tasks 1 and 2.
- Produces: repository navigation with no checked-in documentation HTML or Pages deployment dependency.

- [ ] **Step 1: Remove Pages-only files**

Delete `.github/workflows/pages.yml` and `docs/.nojekyll`.

- [ ] **Step 2: Update project metadata**

Set the `Documentation` URL in `pyproject.toml` to:

```toml
Documentation = "https://github.com/jsonmad/QXYCell/tree/main/docs"
```

- [ ] **Step 3: Update documentation links**

Change the four root README reference links and the focused-guide function-reference links from `.html` to `.md`. Update converted-document links to the new Markdown targets.

- [ ] **Step 4: Validate all local documentation targets**

Run a Python link scan across `README.md` and every `docs/*.md`, ignoring external URLs and fragment-only links, and fail when any relative target does not exist.

- [ ] **Step 5: Verify removal and formatting**

Run:

```bash
test -z "$(find docs -type f -name '*.html' -print)"
test ! -e .github/workflows/pages.yml
test ! -e docs/.nojekyll
! rg -n 'jsonmad\.github\.io/QXYCell|QXYCell_(overview|flowchart|function_reference)\.html|qxy_function_examples\.html|docs/index\.html' README.md docs pyproject.toml .github
git diff --check
conda run -n qxycell ruff check src
```

Expected: no checked-in documentation HTML, Pages configuration, stale HTML links, whitespace errors, or Ruff failures.

- [ ] **Step 6: Commit the implementation**

Run:

```bash
git add -A
git commit -m "docs: replace HTML documentation with Markdown"
```

Expected: one implementation commit containing the Markdown conversion, extracted assets, link updates, and HTML/Pages removals.
