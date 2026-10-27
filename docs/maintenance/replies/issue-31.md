Thank you for including `sessionInfo()`. Your report uses R 4.3.1, Matrix 1.6-5, irlba 2.3.5.1, Signac 1.14.0, Seurat 5.2.1 and SeuratObject 5.0.2. The failure is in Signac's SVD dependency before scMEGA network inference.

This error is consistent with an irlba binary built against an older Matrix C interface. Matrix 1.6-2 changed that interface and requires affected linking packages to be rebuilt. Package version numbers alone do not show which Matrix headers were used to compile irlba. See [Matrix's release notes](https://cran.r-project.org/web/packages/Matrix/news.html) and the [upstream irlba report](https://github.com/bwlewis/irlba/issues/70).

Please restart R and check for multiple library installations:

```r
.libPaths()
packageVersion("Matrix")
packageVersion("irlba")
find.package("Matrix")
find.package("irlba")
```

To preserve the irlba version in your report, rebuild its source in a dedicated library against the Matrix installation you intend to use (a working compiler is required):

```r
lib <- path.expand("~/R/scmega-compat")
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/irlba/irlba_2.3.5.1.tar.gz",
  repos = NULL, type = "source", lib = lib
)
```

Restart R again, put that library first, and test sparse SVD before rerunning `RunSVD()`:

```r
.libPaths(c(path.expand("~/R/scmega-compat"), .libPaths()))
set.seed(31)
a <- Matrix::rsparsematrix(80, 40, density = 0.2)
irlba::irlba(a, nv = 3)$d
```

Verification: on R 4.5.3/Matrix 1.7-5, both the existing irlba 2.3.7 and a separately source-rebuilt irlba 2.3.5.1 passed sparse SVD comparisons against a dense reference and a synthetic Signac TF-IDF → feature selection → SVD pipeline. This verifies the rebuild mechanism on that environment; your original R 4.3.1/Ubuntu binary was unavailable. If it still fails, please share the compilation log and paths above, and confirm the R session was restarted.

No scMEGA algorithm change or blanket Matrix downgrade is indicated. In particular, do not install Matrix 1.7.x into R 4.3.1: that Matrix series requires R >= 4.4.
