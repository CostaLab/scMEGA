# Repository audit

Audit date: 2026-10-27. Baseline: `e6ce4f9` on `main`; clean working tree before work. GitHub API confirms default branch `main`; authenticated viewer has ADMIN access. This audit does not authorize remote writes.

## Architecture and workflow

scMEGA 1.1.0 is an R >= 4.3.0 package with nine R files, two C++ files (Rcpp registration and row-wise correlation), 26 Rd manuals, four bundled gene-annotation RDA files (hg19, hg38, mm9, mm10), four R Markdown vignettes, and a checked-in pkgdown site. No `tests/`, `inst/`, or `data-raw/` existed at baseline. There is no dependency lockfile.

`CoembedData` integrates RNA and gene-activity representations using Seurat anchors; `PairCells` performs greedy or geodesic matching; `CreatePairedObject` transfers paired counts, embeddings and RNA metadata into a multimodal object. Diffusion maps and supervised trajectories feed `GetTrajectory`, TF/gene selection, peak-to-gene correlations and GRN construction. Visualization includes heatmaps and graph layouts. The multiome vignette also uses already paired observations.

Public exports (including the Rcpp wrapper):

`AddTargetAssay`, `AddTrajectory`, `CellPropPlot`, `CoembedData`, `CompareCellProp`, `CorrelationHeatmap`, `CorrelationPlot`, `CreatePairedObject`, `GRNHeatmap`, `GRNPlot`, `GRNSpatialPlot`, `GetCorrelation`, `GetGRN`, `GetPairList`, `GetTFGeneCorrelation`, `GetTrajectory`, `NetCentPlot`, `PairCells`, `PeakToGene`, `PseudotimePlot`, `RunDiffusionMap`, `SelectGenes`, `SelectTFs`, `TopEmbGRN`, `TrajectoryHeatmap`, `TrajectoryPlot`, `rowCorCpp`.

Internal helpers include `getQuantiles`, `centerRollMean`, and `splitAndFetch`. `rowCorCpp` is exported, not internal. Rcpp code checks matrix column counts and upper row-index bounds but deserves separate review of empty/invalid indices. No scientific source was changed during this audit.

## Dependencies and compatibility

DESCRIPTION Imports: Rcpp, ArchR, SeuratObject, dplyr, optmatch, ggplot2, cowplot, ggpubr, destiny, pracma, igraph, glue, ComplexHeatmap, IRanges, Matrix, circlize, ggrepel, matrixStats, nabor, stringr, textshape, tidyr, GenomicRanges, S4Vectors, SummarizedExperiment, ggraph, factoextra. LinkingTo: Rcpp. Suggests: knitr, rmarkdown, testthat >= 3.0.0; testthat edition 3 already configured.

* Seurat is imported in NAMESPACE and used extensively but absent from Imports/Depends. `Remotes: satijalab/seurat@seurat5` is not a dependency declaration. This is an actual package-check ERROR.
* Signac supplies the unqualified `CreateChromatinAssay` call but is absent from DESCRIPTION/NAMESPACE. FNN and data.table are used via `::` but undeclared. Several ArchR/ComplexHeatmap/SummarizedExperiment/genomic accessors are unqualified without corresponding imports; attaching packages interactively may conceal namespace failures.
* Seurat v5 APIs `LayerData` and `JoinLayers` are used. No minimum SeuratObject version is declared; current source cannot be advertised as v4-compatible. Legacy Assay and Assay5 representations must be tested separately; ChromatinAssay must retain genomic ranges/motifs rather than be blindly coerced to Assay5.
* SeuratObject's layer accessor selects cells but Assay5 returns storage order. Blind renaming in CreatePairedObject is the #32 hypothesis. Split layers can return only the first matching counts/data layer (#27); joining or rejecting split layers needs an explicit scientific contract.
* Matrix is a direct dependency; irlba is not called in scMEGA R/C++ source and is a workflow/transitive dependency through Seurat/Signac. SingleCellExperiment is neither declared nor directly used in package source; its presence in user sessions does not make it a direct dependency. GenomicRanges and SummarizedExperiment are genuine direct dependencies.
* #31 reports a Matrix/irlba compiled-interface mismatch. Rebuild binary-linked packages against a coherent Matrix installation rather than rewriting analysis algorithms. No dependency upgrades were performed.
* Missing version constraints should be established from validated combinations, not guessed. The README requests historical `seurat5` branches for Seurat and Signac; the install vignette differs. ArchR also needs explicit source provenance. Recommend a separate maintainer-reviewed metadata/install change.
* Deprecated igraph aliases (`shortest.paths`, `clusters`) and direct slot access warrant compatibility tests. Do not replace them solely to hide warnings.

## Scientific risks identified by inspection (not confirmed defects)

Cell/metadata/reduction ordering, singleton matrix dimension dropping in GetTrajectory, empty trajectory bins and normalization, motif-symbol correspondence, gene annotation ordering in PeakToGene, and absent/duplicate TF identifiers merit deterministic invariants. `RunDiffusionMap` appears to ignore non-NULL `k` and uses 30; investigate separately. No conclusions about biological validity of these paths have been established by this audit.

## Baseline validation

Local host: macOS arm64, R 4.5.3, Seurat 5.5.0, SeuratObject 5.4.0, Matrix 1.7-5, irlba 2.3.7. Signac, testthat, ArchR, GenomicRanges and SummarizedExperiment were absent.

* `R CMD INSTALL --library=/tmp/scmega-baseline/lib .`: FAILED, 13 Imports unavailable: ArchR, optmatch, ggpubr, destiny, pracma, ComplexHeatmap, circlize, nabor, textshape, GenomicRanges, SummarizedExperiment, ggraph, factoextra.
* `R CMD build --no-build-vignettes --no-manual .`: PASSED, produced scMEGA_1.1.0.tar.gz (kept outside repository).
* `R CMD check --no-manual --no-build-vignettes --output=/tmp/scmega-baseline scMEGA_1.1.0.tar.gz`: FAILED, 1 ERROR at dependency checks, also missing suggested testthat and undeclared Seurat. Initial sandbox run also could not access package indexes; network retry recorded separately.
* There were no existing tests to run. This is not a passing baseline. Full vignette rendering was not attempted because biological datasets and dependencies are missing.

## CI and recommendations

Only `.github/workflows/render-rmarkdown.yaml` existed. It runs on every push, installs moving GitHub dependencies and the remote scMEGA instead of the checked-out candidate, builds the website, then commits and pushes. It has no explicit permissions boundary or package checks. Automatic publication branches must not trigger this push-back path; this task will contain it with a manual trigger and read-only permissions.

Add a minimal testthat runner and synthetic pairing fixtures, a separate R package check, and isolated maintenance security tests. Start with Linux R 4.5 and v5-oriented dependencies as a *candidate*, not a supported matrix. Keep CI red on real metadata/dependency problems. Expand to R minimum/current and compatible Bioconductor releases only after installation and baseline checks work. Preserve README, algorithms, public defaults, genome conventions and bundled reference data.

References: [SeuratObject layers](https://satijalab.github.io/seurat-object/reference/Layers.html), [r-lib actions](https://github.com/r-lib/actions), [Matrix release notes](https://cran.r-project.org/web/packages/Matrix/news.html). Official layer docs specify cells to include, not guaranteed requested ordering; local implementation inspection and reproduction are decisive for #32.
