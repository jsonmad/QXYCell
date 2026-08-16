# Staged workflow example rules

- Keep one QXYCell stage per numbered script.
- Read all user-editable settings from `config.py`.
- Keep the two threshold-source scripts separate and explicit; never add fallback selection.
- Use generic placeholder paths and context only. Do not add private data or research paths.
- Keep every script import-safe with a `main()` guard.

