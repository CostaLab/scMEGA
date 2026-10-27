# fix: preserve trajectory matrix dimensions (#25)

Related: #25. Local draft; no remote PR.

GetTrajectory drops singleton cell/feature dimensions before rowMeans, and single-bin aggregation can return a vector. GetCorrelation likewise fails with one shared TF. Preserve those matrix dimensions and give an explicit error for empty trajectory bins instead of undefined means.

The public API, nonempty-bin means/normalization and correlation/FDR formulas are unchanged. Empty-bin behavior is intentionally tightened and documented; callers must choose a justified grouping rather than fill missing means artificially. Long smoothing windows, zero-time endpoint policy and other input edge cases are outside this correction.

Regression tests fail on original source and pass after correction: 29 expectations, no warnings/failures/skips. Installed-package tests also pass in a scratch copy with separately proposed metadata commit ec24fdf. Installation/build pass; R CMD check completes with six existing warnings and three notes, so CI is not clean. Original reporter data were unavailable; this confirms matching reproducible defects, not the exact failing bin in that dataset.

Branch codex/issue-25-trajectory, commit 98ff9ad, based on main e6ce4f9. The metadata proposal is a separate review item. No merge or publication is authorized by this draft.
