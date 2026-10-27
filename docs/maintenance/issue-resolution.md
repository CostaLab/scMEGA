# Phase II issue resolution

Started 2026-10-27; continued 2026-10-28. Current GitHub list refreshed: 12 open issues. Existing commits a9bea1e, 6c3ff52 and dde1d7c are preserved. Main remains e6ce4f9. No remote writes are authorized or performed. Code fixes use independent worktrees based on origin/main; this reporting branch is based on the completed Phase I branch.

Investigations proceed in the requested order: #25, #31, #24, #27, #33, #23, #21, #18, #17, #30, #29. #32 is rechecked for preservation. Dependency installation uses the isolated library and targeted repairs of observed blockers; it did not delay the sequential source-grounded issue investigations.

| Issue | Category | Priority | Status | Resolution | Tests | Branch / Reply |
|---|---|---|---|---|---|---|
| #32 | Confirmed bug | P0 | Existing local fix preserved and reverified | Named RNA/ATAC pairing | 117 source and installed expectations pass (with metadata); 2 warnings | codex/issue-32-createpairedobject; 6c3ff52; replies/issue-32.md |
| #25 | Confirmed bug | P1 | Fixed and source-tested locally | Preserve matrix dimensions; reject empty bins explicitly | 29 source and installed expectations pass (with metadata); original edge cases fail | codex/issue-25-trajectory; replies/issue-25.md |
| #31 | Compatibility issue | P1 | Diagnosed; compatibility answer drafted | Rebuild irlba against active Matrix | Current and rebuilt 2.3.5.1 SVD/Signac tests pass | replies/issue-31.md; no code branch |
| #24 | Compatibility issue | P1 | Partial compatibility verified; answer drafted; needs reproduction | Distinguish missing/split layers, assay classes, dimensions and imports | Layer/count/data/ChromatinAssay probes pass | replies/issue-24.md; metadata proposal separate |
| #27 | Compatibility issue | P1 | Verified preprocessing solution plus local compatibility fix | Explicit RNA layer preparation; preserve legacy feature selection | 17 source and installed expectations pass (with metadata); real 220-cell pipeline | codex/issue-27-coembed; replies/issue-27.md |
| #33 | Usage question | P1 | Diagnosed; answer drafted; mapping policy/data needed | Exact labels are not species-aware identifiers | Signac metadata and synthetic matching verified | replies/issue-33.md; no speculative remapping |
| #23 | Compatibility issue | P1 | Layer diagnosis verified; answer drafted; original object needed | Normalize intended RNA input; preserve peak motif class; activity can be ordinary Assay | Missing-data and corrected trajectory probes pass | replies/issue-23.md; no code branch |
| #21 | Usage question | P2 | Verified answer drafted | Weak component plus highlights; correct vertex-name/type assignment | Error reproduced; plot built/rendered and label inspected | replies/issue-21.md; no code branch |
| #18 | Feature request | P1 | Evidence verified; answer drafted; policy decision needed | Curated species-aware mapping across the full pipeline | Reuses verified #33 motif diagnostics | replies/issue-18.md; no code branch |
| #17 | Feature request | P3 | Scope verified; answer drafted; annotation/API decision needed | Four bundled genomes only; custom annotations require validated extension | Unsupported rn6 rejected; 3 mm10 links/C++ correlations verified | replies/issue-17.md; no code branch |
| #30 | Usage question | P3 | Source-grounded scientific answer drafted | Trajectory-conditioned scope; no skip-trajectory switch | Synthetic GetGRN evidence joins/support counts verified | replies/issue-30.md; no code branch |
| #29 | Documentation issue | P3 | Existing diagnosis preserved; answer drafted | Correct current arguments and original-count inputs | Screenshot/argument error already verified; paired suite rechecked | replies/issue-29.md; existing dde1d7c |

## #25 — trajectory dimensions (processed first)

The report and maintainer comment were refreshed and read. A matching error was reproduced in the real GetTrajectory function with sparse Seurat inputs. Singleton cells or a single feature collapse `data.use[, cell_names]` to a vector. A one-bin Reduce(cbind, ...) also returns a vector. GetCorrelation similarly drops dimensions for a single shared TF. Empty bins have undefined means and previously failed later or propagated NaN when normalization was disabled.

Independent worktree `/tmp/scmega-issue-25`, based on main e6ce4f9. Changes: `R/trajectory_analysis.R`, `man/GetTrajectory.Rd`, `tests/testthat.R`, `tests/testthat/test-trajectory-dimensions.R`. No dependency on #32. Added drop=FALSE, matrix-preserving column binding and an explicit empty-bin error; no public defaults/statistical method changed. Existing nonempty multi-cell results are checked against explicit means and normalization. The invalid metadata and no-common-TF errors already existed and were preserved.

Command: `R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/run-issue25-tests.R /tmp/scmega-issue-25`. Original: four test errors; ordinary means and existing validation pass. Corrected: 29 expectations, zero failures/errors/skips/warnings. Source harness explicitly binds the intended SeuratObject Assays accessor because attaching SummarizedExperiment masks that name; no constructor/accessor is mocked. Installed-package validation is separate. Logs: phase2-evidence/issue-25-before.log and issue-25-after.log.

Remaining: reporter object absent, so the precise failing bin/feature is unknown. Smoothing windows longer than the number of bins, zero-time endpoint exclusion, non-dividing groupEvery values and absent assay cells need explicit analysis contracts; no speculative interpolation, trajectory assignment, endpoint change or motif remapping was introduced. Empty-bin input now fails deliberately instead of producing undefined means, which is documented. Affected previously failing cases can be rerun after review; valid nonempty outputs remain numerically unchanged. Maintainer review required; reply drafted, nothing published.

#25 local commit: `98ff9ad5284ce9b4fe285baf8759ba07bb9d04e7`. Actual package validation: build exit 0; install/check exit 1, missing dependencies and undeclared Seurat. No remote PR.

## #31 — Matrix / irlba (processed second)

Full report: R 4.3.1, Matrix 1.6-5, irlba 2.3.5.1, Signac 1.14.0, Seurat 5.2.1, SeuratObject 5.0.2. No comments or attachments. scMEGA source does not directly call irlba; the failing RunSVD is Signac's function. The missing C callable is consistent with a stale compiled interface, supported by Matrix's official 1.6-2 ABI/rebuild notes. The original binary's compilation history is unavailable, so the exact library mismatch remains inferred.

Downloaded the official CRAN irlba 2.3.5.1 source archive and rebuilt it into /tmp/scmega-irlba2351-lib without replacing the existing irlba. `R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-31-matrix-irlba.R` and the same command with `R_LIBS_USER=/tmp/scmega-irlba2351-lib:/tmp/scmega-r-lib` both pass: singular values match dense SVD within 1e-6; relative residual 5.523446e-16; Signac TF-IDF/feature selection/RunSVD returns finite 60-by-5 embeddings with correct cell IDs. R 4.5.3, Matrix 1.7-5; tested irlba versions 2.3.7 and rebuilt 2.3.5.1. Exact R 4.3.1/Ubuntu environment not reproduced. No source change or issue branch needed. Reply gives version-preserving source-rebuild instructions and asks for logs if the error persists; no blanket downgrade/upgrade recommendation. Not published.

## #24 — Seurat v5 (processed third)

The refreshed issue only says GetTrajectory fails with new Seurat; no comments, versions, traceback or attachment. Current source is already migrated to LayerData, but this alone does not establish compatibility. Runtime namespace audit finds LayerData and JoinLayers are not exported by imported Seurat 5.5.0 and not explicitly imported from SeuratObject in scMEGA. That metadata defect has a separate proposal; do not mix it into an algorithm rewrite.

`R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-24-seurat-layers.R` passes: counts and normalization preserve expected values under legacy Assay and Assay5; split counts return only the first matching layer (warning asserted), explicit JoinLayers recovers all cell-named counts; counts-only objects lack data; ChromatinAssay layer access preserves ranges and class. GetTrajectory's dimension correction belongs to #25, pairing to #32, split transfer reference investigation to #27. SelectTFs/GetTFGeneCorrelation require ChromatinAssay motif metadata; do not cast it away. No broader code fix without an exact failing case. Reply drafted requesting reproducible inputs, not published. No new issue branch needed.

## #27 — split reference data and legacy merging (processed fourth)

Refreshed full report/sessionInfo read, no comments/attachments. Reported 2,976 data columns vs 4,224 anchor reference cells plus multiple counts/data layer warnings. Real Seurat pipeline on synthetic 100-gene, 120-RNA/100-ATAC inputs reproduced the same error (60 vs 120). Joining intended disjoint RNA layers and performing ordinary RNA normalization/variable-feature selection completed anchor transfer/PCA/UMAP on all 220 cells. The first fixture omitted variable-feature selection and failed later; that was an input-preprocessing omission, not evidence to change the algorithm. The corrected fixture explicitly supplies the reference selection. Original UMAP weighting also requires checking dimensions: a two-dimensional embedding cannot supply default dims=1:30.

The legacy Assay comparison revealed a separate reproducible compatibility failure: unconditional JoinLayers has no Assay method. Guarding that exposed legacy merge dropping the user-provided variable features. The minimal branch checks the class and restores the supplied reference variable-feature set only when the merged object has none. It never invents a new variable-feature selection. Split/absent normalized data layers now fail before anchor calculation with preparation instructions rather than selecting a first batch silently. No automatic joining of ambiguous input layers or changes to normalization, PCA/UMAP defaults, feature-selection method or public API.

Independent worktree /tmp/scmega-issue-27 based on main e6ce4f9; no dependency on #25/#32. Files: R/pair_cells.R, man/CoembedData.Rd, tests/testthat.R, tests/testthat/test-coembed-layers.R. `R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/run-issue27-tests.R /tmp/scmega-issue-27`: original has two test errors; final 17 expectations, zero failures/errors/skips, five upstream warnings (missing imputed RNA counts, large PCA rank fraction, UMAP implementation notice). Observed normalized RNA data, ATAC counts, identifiers, reference variable-feature selection and finite embeddings are checked. Full diagnostic examples/logs retained, including failed intermediate attempts. Not a validation of plant reference genomes or biology. Maintainer review required; reply drafted, nothing published.

Local 27 commit: `7fbcce628ee14884dfcf4e588a6bfcc70bc3d2a6`. Build passed; install/check remain blocked by ArchR (the metadata proposal removes the undeclared-Seurat dependency error).

Local metadata commit: `ec24fdf3943128b37fb99bc51148c3eb86f6b6bc`. Build passed; install/check remain blocked by ArchR (the metadata proposal removes the undeclared-Seurat dependency error).

## #33 — mm10 downstream analysis (processed fifth)

Full issue read; no comments/attachments. The shown AddMotifs result is assigned to seurat_LTHSC while SelectTFs receives seurat_obj, a separate possible usage mistake. SelectTFs has no genome argument and matches motif labels to RNA names exactly via GetCorrelation; the mm10 sequence input to AddMotifs does not implement gene-ID/orthology mapping. Case differences, composite motifs, duplicate TF labels and RNA Ensembl-vs-symbol identifiers require explicit scientific mapping. Species/build and chromosome naming still matter upstream and in PeakToGene, but are not consulted by this exact-string matching step. SelectGenes defaults to hg38, so mm10 must be explicitly selected downstream.

`R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-33-motif-identifiers.R` verifies real Signac metadata, retained Spi1/Actb vs lost Gata1, duplicate/composite ambiguity under title casing, and why motif IDs must be checked before positional label assignment. Bundled mm10 gene symbols and seqnames inspected. These are component/identifier tests, not full SelectTFs, genome scanning or mouse biological validation (this component-stage test preceded ArchR installation; reporter data remain unavailable). The row-order risk is recorded for follow-up, not claimed as the original report's demonstrated cause. No code branch: automatic capitalization, duplicate removal or splitting complexes would change biology without an agreed policy. Reply requests curated examples/provenance and consistent object use. Maintainer mapping/API decision required; not published.

## #23 — missing layers and assay conversion (processed sixth)

Full report and all five returned comments read. The old GetAssayData traceback is not the current source and does not establish a rat genome error. Current SelectTFs requests data layers in both RNA and TF activity; class conversion cannot create those data. The comment's SoupXRNA adaptation is already supported by rna.assay. Peak ChromatinAssay must retain motifs; chromVAR activity need not itself be a ChromatinAssay. return.heatmap=FALSE still computes an internal TrajectoryHeatmap, so it is not an ArchR-free workaround.

`R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-23-missing-layers.R` passes diagnostic assertions: counts-only RNA fails (recorded missing-layer warning), intended NormalizeData restores full GetTrajectory output without casting, and a data-only ordinary Assay retains signed activity means. Uses nonempty multi-cell bins to avoid conflating #25. Full SelectTFs/rat data and the reporter's old versions are unavailable; request sessionInfo, classes/layers and current body/traceback. No new code branch or arbitrary downgrade; reply drafted, not published.

## #21 — network connectivity and labels (processed seventh)

Read all four comments and inspected the attached network image. NetCentPlot's focus layout requires connectivity; the reporter already confirmed component selection works. Their follow-up type assignment compares numeric vertex positions with a string TF name and uses TF rather than the exact TF/Gene label required by the default source branch. The existing highlights argument provides explicit labels independently of that centrality/type filter.

`R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-21-network-labels.R` reproduces the disconnected-graph error, selects the correct weak component by vertex name, builds/renders the plot and verifies RUNX1 is labeled. The original five-node graph remains unchanged; displayed component has three nodes. Rendered PNG visually inspected. No source change or new branch. Reply explains display scope vs biological inference; not published.

## #18 — mouse motif matching (processed eighth)

Read the report and six comments. The maintainer explicitly noted human-only testing in that discussion; the proposed title-case workaround is not a general orthology/alias map. Source inspection confirms the same positional label assignment in SelectTFs and GetTFGeneCorrelation; GetGRN requires matching TF names in motif.matching columns and correlation rows. Species mapping must preserve stable motif IDs and define duplicate/complex handling consistently, including the hypothesis set for FDR. The verified #33 diagnostic is reused instead of repeated installation or new speculative code. No mapping API is designed without maintainer review and example provenance. #18/#33 are related but not closed as duplicates; #33 adds composite motifs/object-use concerns. Reply drafted, not published.

## #17 — custom genomes and later TF indexing (processed ninth)

Refreshed issue/two comments match the previously read full content. Current PeakToGene accepts four literal assembly names; SelectGenes forwards genome and cannot make a local GRanges annotation effective by itself. The macaque replacement shown in the comment uses placeholder correlation values/random heatmap rather than real peak-to-gene calculations. A proper extension requires assembly/coordinate/TSS/gene-identity conventions, duplicate handling and annotation provenance, with no nearest-species substitution or arbitrary chromosome prefix changes across builds.

`R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-17-genome-contract.R`: rn6 rejected as expected; real bundled mm10 annotation, GenomicRanges interval operations and compiled rowCorCpp produce three links whose correlations match direct named-row calculations within 1e-10 after gene-row permutation. No real dataset download or algorithm changes. SelectGenes hg38 default and PeakToGene hg19 default verified. Full custom-species support is a feature decision, not a string addition.

The secondary TF indexing error was inspected separately in GetTFGeneCorrelation/TrajectoryHeatmap: requested names may differ after one-sided remapping or be absent after filtering; reporter-modified function/object unavailable. Reply asks for requested IDs and post-filter activity row names rather than silently dropping TFs. No speculative fix or code branch; maintainer review needed for a new annotation/mapping interface. Not published.

## #30 — trajectory scope and scientific interpretation (processed tenth)

Full report/two comments refreshed and reread. GetTrajectory excludes NA metadata and groups cells into intervals before normalization/smoothing. SelectTFs, SelectGenes and GetTFGeneCorrelation always use this path; trajectory.name only selects a metadata column. No trajectory-free switch exists. A CD4 trajectory cannot establish networks for every PBMC lineage; arbitrary cell-order pseudotime would change the statistical question.

`R_LIBS_USER=/tmp/scmega-r-lib Rscript docs/maintenance/examples/issue-30-grn-meaning.R` verifies actual GetGRN joins: two peak rows support TF1/G1, one unmatched TF2/G1 row retains NA correlation, and a supplied negative TF2/G2 correlation remains negative. The function does not itself select significant edges or validate causality. n_peaks is a support-row count, not a regulatory probability; input preparation and duplicate mapping matter. No biological all-PBMC analysis or statistical-method change. Reply distinguishes predictions from experimental regulation and requests a separate validated design for trajectory-free inference. No code branch; not published.

## #29 — existing documentation diagnosis (processed eleventh)

Current body/attachment URL unchanged from Phase I. Reused the inspected immutable screenshot and existing unused-arguments reproduction rather than duplicating the completed investigation. Source signature and corrected docs remain intact on the reporting base; the 117 pairing expectations were rerun this phase. Reply provides obj.coembed/obj.rna/obj.atac/rna.assay/atac.assay, explains df.pair and the distinct original count objects, and states the local/unpublished status. No new branch, source change or duplicate fix. Original dde1d7c preserved.

## #32 — preservation check

Original a9bea1e/6c3ff52/dde1d7c commits preserved. Reporting branch descends from them; independent #25/#27/metadata branches start from upstream main and deliberately do not contain unrelated fixes. `R_LIBS_USER=/tmp/scmega-r-lib SCMEGA_REQUIRE_INTEGRATION=true Rscript scripts/maintenance/run_pairing_tests.R .` passes all 117 expectations, zero failures/errors/skips and two expected dense-input conversion warnings. The 22 security tests also pass and automation files have not been changed. Reply drafted. No branch was pushed, merged or published.

## Final review and publication status

All 12 open issues have an outcome and a reply draft. The final GitHub refresh confirms the same issue numbers and update timestamps; no new report/comment was missed during this session. All replies remain local. No branch push, remote PR/comment/label, issue closure, branch merge or automation enablement occurred.

Independent code commits, each based on upstream main e6ce4f9:

| Branch | Commit | Files | Review scope |
|---|---|---|---|
| codex/issue-25-trajectory | 98ff9ad | R/trajectory_analysis.R; man/GetTrajectory.Rd; tests/testthat.R; tests/testthat/test-trajectory-dimensions.R | Dimensions and explicit empty-bin handling; no statistical-method change |
| codex/issue-27-coembed | 7fbcce6 | R/pair_cells.R; man/CoembedData.Rd; tests/testthat.R; tests/testthat/test-coembed-layers.R | RNA layer contract, legacy class handling, retention of supplied reference features |
| codex/dependency-metadata | ec24fdf | DESCRIPTION; NAMESPACE; R/imports.R | Separate runtime dependency/import review; generated namespace, no algorithms |

The reporting branch codex/phase2-resolution descends from dde1d7c and preserves #32. Its new files are reports, reply drafts, executable verification examples and evidence. The independent worktrees are under /tmp; their committed Git branches persist even if those temporary directories are removed. The fixes have not been combined: maintainer integration testing is required, particularly for #27/#32 touching different functions in R/pair_cells.R and for the metadata proposal enabling package namespace execution.

#25/#27/#32 are suitable for human review as locally verified corrections, with the recorded validation limits. Metadata is a separate review proposal. No fix is approved for automatic publication or merge. #24/#23 still need exact current failing objects/tracebacks; #33/#18 need curated mappings and a duplicate/composite policy; #17 needs an annotation interface/assembly provenance and information for its secondary TF error. #21/#30/#29 have verified/source-grounded usage explanations ready for review as replies. #31's source-rebuild mechanism is tested, but the original Ubuntu/R 4.3.1 binary was unavailable.

A complete green installed-package check, a coherent portable dependency stack, and hosted automation/security/model/PR validation remain separate acceptance items. No existing security check was weakened. Logs distinguish source tests, actual package checks and caught diagnostic errors. No untested biological result or original dataset has been claimed validated.

## Final installed-package verification (2026-10-28)

The environment was subsequently completed through isolated system-library and compiler-standard repairs; see dependency-audit.md. All declared direct Imports load. Original reporting-base package installation/build now pass; its check still fails on the undeclared Seurat dependency. The metadata-only proposal completes check with 6 WARNINGs/3 NOTEs, no errors. Each separate fix plus metadata was tested in a disposable filesystem copy (not a Git merge): installed #25 tests PASS 29/WARN 0, #27 PASS 17/WARN 5, #32 PASS 117/WARN 2, all FAIL 0/SKIP 0. Their built-tarball checks complete with the same 6 warnings/3 notes. This supersedes earlier missing-dependency blockers; it does not turn a warning-bearing check into clean CI success. Source and installed test modes are explicitly distinguished in the retained logs.

Remaining warning/note categories are already present in the metadata-only baseline. Hosted workflow/security/model/PR behavior remains untested and automatic publication stays disabled. The local code fixes and metadata proposal are ready for human draft-PR review, with warning cleanup/integration testing still required before merge/rollout. Replies #24/#23 and mapping/custom-genome requests remain data/design dependent, not declared resolved simply because a reply exists.

Full commands/results: [Phase II validation](phase2-validation.md). All responses: [reply index](replies/README.md). Protected-file proposal and environment details: [dependency audit](dependency-audit.md). Local review descriptions: [#25](pr-drafts/issue-25.md), [#27](pr-drafts/issue-27.md), [metadata](pr-drafts/dependency-metadata.md).

The new local issue/metadata branches have no upstream tracking, preventing a plain upstream push from targeting origin/main. No push was performed.

Complete path/commit inventory: [Phase II changed files](phase2-changed-files.md).
