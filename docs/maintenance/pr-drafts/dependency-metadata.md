# chore: declare runtime dependencies and explicit namespace imports

Local proposal requiring separate maintainer review of protected DESCRIPTION/NAMESPACE changes. No remote PR.

Seurat is imported but undeclared; Signac, FNN, data.table and viridis are used at runtime but undeclared. Several runtime helpers, including LayerData/JoinLayers, are not available through the existing imports. Declare the actual providers, specify v5 API minima, and replace the historical Seurat development remote with the reviewed ArchR source commit. Add explicit roxygen namespace directives and regenerate NAMESPACE with recorded roxygen2 8.1.1.

No scientific algorithm, function signature/default, package version or existing site-library version is changed. Version minima describe used APIs, not a tested compatibility matrix. Vignette-only dependencies/data provisioning still need separate decisions.

After isolated environment repair, installation/build pass and the metadata-only package check completes with zero errors, six warnings and three notes. Scratch copies with each independent #25/#27/#32 correction pass installed tests (29/17/117 expectations). Warning/note categories remain and are not suppressed. Private ArchR calls and the mixed local Bioconductor stack remain portability concerns; hosted CI is unverified.

Branch codex/dependency-metadata, commit ec24fdf, based on main e6ce4f9. Files: DESCRIPTION, NAMESPACE, R/imports.R. This proposal cannot use the automatic patch publication path.
