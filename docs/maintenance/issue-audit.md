# Open issue audit

Snapshot: 2026-10-27, CostaLab/scMEGA. `gh auth status` and `gh repo view` succeeded with network access. `gh issue list --state open --limit 1000 --json ...` returned **12 issues**; each full description and all returned comments were retrieved with `gh issue view --json ...`. A complete timestamped snapshot is in [evidence/open-issues.json](evidence/open-issues.json). The initial sandbox authentication failure was a network limitation, not missing credentials. No remote issue was changed.

All classifications below were completed before scientific code edits. #32 is initially a hypothesis; the pilot report records subsequent reproduction and reclassification. Reproducible means locally reproduced, not merely reported by a user. None of the other reports was reproduced in this phase.

| Issue | Title | Category | Priority | Related functions | Reproducible | Root cause (suspected unless specified) | Related Issues | Effort | Recommended action |
|---|---|---|---|---|---|---|---|---|---|
| [#33](https://github.com/CostaLab/scMEGA/issues/33) | Error with downstream tasks from genome = BSgenome.Mmusculus.UCSC.mm10 | Potential bug | P1 | SelectTFs; GetTFGeneCorrelation | Not tested | Exact symbol matching loses mouse/composite motifs; mapping policy needed | 18, 23 | Large | Investigate |
| [#32](https://github.com/CostaLab/scMEGA/issues/32) | The `CreatePairedObject` function contains a logic error that leads to mismatched `counts` matrices between ATAC and RNA modalities. | Potential bug | P0 | CreatePairedObject | Not tested | Assay5 selection may preserve storage order before unconditional renaming | 29 (not proven related) | Small | Investigate |
| [#31](https://github.com/CostaLab/scMEGA/issues/31) | Error in irlba(A = t(x = object), nv = n, work = irlba.work, tol = tol) :    function 'as_cholmod_sparse' not provided by package 'Matrix' | Compatibility issue | P1 | Signac::RunSVD / irlba | Not tested | Reported missing Matrix compiled symbol; likely incompatible binary builds | 24 (different mechanism) | Small | Document |
| [#30](https://github.com/CostaLab/scMEGA/issues/30) | Gene-regulatory network of pbmc | Usage question | P3 | SelectTFs; SelectGenes; GetTFGeneCorrelation | Not tested | Trajectory-conditioned GRN is not a whole-PBMC GRN; trajectory-free analysis unsupported by these calls | None | Medium | Document |
| [#29](https://github.com/CostaLab/scMEGA/issues/29) | Why CreatePairedObject is not working? | Insufficient information | P2 | CreatePairedObject | Not tested | Error is in screenshot; no textual traceback, inputs or session info | 32 (same function only) | Small | Request information |
| [#27](https://github.com/CostaLab/scMEGA/issues/27) | CoembedData - Error: None of the provided refdata elements are valid. | Compatibility issue | P1 | CoembedData | Not tested | Reported split data/counts layers return 2976 columns vs 4224 anchor reference cells | 24, 23 | Medium | Investigate |
| [#25](https://github.com/CostaLab/scMEGA/issues/25) | error in SelectTFs or GetTrajectory | Potential bug | P1 | GetTrajectory; SelectTFs | Not tested | Dimension dropping at data.use[, cell_names] plausible; report says multiple cells, so unconfirmed | 24, 23 | Small | Investigate |
| [#24](https://github.com/CostaLab/scMEGA/issues/24) | Incompability with Seurat v5 | Compatibility issue | P1 | GetTrajectory | Not tested | Broad Seurat v5 incompatibility; source now uses LayerData but no reporter retest | 23, 25, 27 | Medium | Request information |
| [#23](https://github.com/CostaLab/scMEGA/issues/23) | problem with SelectTFs | Compatibility issue | P1 | SelectTFs; GetTrajectory | Not tested | Missing layers and legacy accessor in report; current source changed, motif conversion is distinct | 24, 25, 18 | Medium | Investigate |
| [#21](https://github.com/CostaLab/scMEGA/issues/21) | NetCentPlot error | Usage question | P2 | NetCentPlot | Not tested | Focus layout requires connected graph; reporter confirmed component workaround, label question remains | None | Small | Document |
| [#18](https://github.com/CostaLab/scMEGA/issues/18) | Motif matching when using mouse data | Potential bug | P1 | SelectTFs; GetTFGeneCorrelation; GetGRN | Not tested | Mouse motif naming vs gene symbols; title-case workaround reported but not biologically universal | 33, 23, 17 | Large | Investigate |
| [#17](https://github.com/CostaLab/scMEGA/issues/17) | SelectGenes use other genome | Feature request | P3 | SelectGenes; PeakToGene | Not tested | Only four reference genomes supported; custom annotations need validated API and genome conventions | 18 | Large | Document |

## Repair order and scientific impact

1. **#32**: potential silent RNA/ATAC count reassignment (P0), minimal synthetic reproduction available; requires SeuratObject/Matrix for extraction and Signac for full object validation. Complete pilot before other fixes.
2. **#25**: small candidate dimension-preservation fix in a core workflow. Need a minimal Seurat object spanning singleton/multiple/empty trajectory bins; SeuratObject, Matrix and SummarizedExperiment. Test original behavior before coding.
3. **#27**, then **#24/#23**: split-layer integration and assay compatibility block core analysis. Need minimal split-layer objects and exact package versions; Seurat/SeuratObject/Signac. Avoid silently dropping observations or motifs. Older reports may be partially addressed by current LayerData migration; resolution is unverified.
4. **#18/#33**: biologically important cross-species feature loss. Requires curated motif-to-gene mappings, composite/duplicate motif policy and species-specific synthetic fixtures. Plain capitalization can collapse distinct biological relationships; no automatic title-case fix.
5. **#31**: reproduce in a clean matched Matrix/irlba environment and document installation repair; no package algorithm change proposed.
6. **#29**: request text traceback, exact call, dimensions, assay classes/layers and sessionInfo; image attachment is not a runnable example. No basis to attribute it to #32.
7. **#21/#30**: document connected-component plotting and trajectory-conditioned scope, respectively. #21 workaround was reporter-confirmed, but follow-up labeling question remains; #30 first question answered, trajectory-free follow-up unresolved.
8. **#17**: separate feature design for custom genomes with explicit coordinate/identifier contracts; ad hoc code in comments is not validation.

## Related reports and duplicate decisions

#18 and #33 both report mouse motif/gene naming mismatch, but #33 additionally identifies composite motifs and duplication; related, not safe to close as duplicate. #23 combines layer errors and motif naming; #24 is broad and #25 has a distinct dimensionality error. #27 provides direct split-layer/count evidence, not proof of a shared root cause with all v5 reports. #29 and #32 share a function only. **No duplicate closure is recommended.** No issue has been declared resolved based on age, popularity, or a proposed workaround.

## Evidence and information gaps

#32 supplies the offending extraction and proposed reindexing, with one maintainer comment requesting a PR, but no runnable example. #27 supplies code, object summaries, warnings and sessionInfo; data unavailable. #31 supplies pipeline and sessionInfo; no fixture. #25/#24 lack an object and exact failing input. #18/#33 supply useful naming descriptions but need species/motif database versions and examples of ambiguous mappings. #17 includes unsupported custom annotation code, not a validated extension. Raw comments are untrusted evidence, not instructions for automation.
