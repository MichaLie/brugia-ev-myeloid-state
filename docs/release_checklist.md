# Release Checklist

Use this checklist before citing a repository version or archiving a public release.

## Local verification

- Confirm the working tree contains only intentional package changes.
- Rebuild selected perturbation outputs where the local dependencies are available.
- Rebuild derived summary tables with `scripts/build_summary_tables.py`.
- Rebuild claim-support summaries with `scripts/build_claim_support_summary.py`.
- Rebuild figure assets with `scripts/build_figures.py`.
- Update `environment/software_versions_verified.txt` if verification is repeated in a different environment.
- Confirm `docs/public_data_manifest.md`, `docs/derived_table_manifest.md`, and `docs/module_definitions.md` match the analysis methods.
- Confirm no drafts, publication-package files, literature PDFs, or non-public data are included.

## Citation metadata

- Record the exact Git commit hash cited with the analysis package.
- Confirm `LICENSE.md` and `CITATION.cff` remain current for the release.
- Add or update repository citation metadata after the final author list, release version, and archive DOI are confirmed.
- Create a GitHub release for the submitted repository state.
- Archive the release with a persistent DOI provider such as Zenodo or Figshare.
- Update the associated data-availability statement with the release DOI when available.

## Post-release guardrails

- Do not rewrite the submitted release history.
- Put later improvements in a new commit and, if needed, a new release.
- Keep evidence-bounded language: selected analyses are rerunnable from cached public inputs, while other layers are fixed derived evidence tables for claim inspection and figure regeneration.
