trajectory_fixture <- function(nfeatures = 2L, times = c(10, 20, 30, 40, 60, 70, 80, 90), legacy = FALSE) {
  counts <- Matrix::Matrix(matrix(seq_len(nfeatures * length(times)), nfeatures,
    dimnames = list(paste0("gene", seq_len(nfeatures)), paste0("cell", seq_along(times)))), sparse = TRUE)
  make <- if (legacy) SeuratObject::CreateAssayObject else SeuratObject::CreateAssay5Object
  object <- SeuratObject::CreateSeuratObject(make(counts = counts), assay = "RNA")
  object$Trajectory <- times
  list(object = object, counts = counts)
}

test_that("trajectory means retain cells, genes and values across assay classes", {
  for (legacy in c(TRUE, FALSE)) {
    f <- trajectory_fixture(legacy = legacy)
    result <- GetTrajectory(f$object, assay = "RNA", groupEvery = 25,
      smoothWindow = NULL, log2Norm = FALSE, scaleTo = NULL)
    expected <- vapply(list(1:2, 3:4, 5:6, 7:8),
      function(i) Matrix::rowMeans(f$counts[, i, drop = FALSE]), numeric(2))
    dimnames(expected) <- list(rownames(f$counts), paste0("T.", c(0,25,50,75), "_", c(25,50,75,100)))
    expect_equal(SummarizedExperiment::assay(result), expected)
    expect_equal(SeuratObject::LayerData(f$object, assay = "RNA", layer = "counts"), f$counts)
    normalized <- GetTrajectory(f$object, assay = "RNA", groupEvery = 25, smoothWindow = NULL)
    expect_equal(SummarizedExperiment::assay(normalized), log2(sweep(expected, 2, colSums(expected), "/") * 10000 + 1))
  }
})

test_that("singleton cells and features remain two dimensional", {
  for (legacy in c(TRUE, FALSE)) {
    f <- trajectory_fixture(times = c(10,30,60,90), legacy = legacy)
    result <- GetTrajectory(f$object, assay = "RNA", groupEvery = 25,
      smoothWindow = NULL, log2Norm = FALSE, scaleTo = NULL)
    expect_equal(unname(SummarizedExperiment::assay(result)), unname(as.matrix(f$counts)))
    f <- trajectory_fixture(nfeatures = 1L, legacy = legacy)
    result <- GetTrajectory(f$object, assay = "RNA", groupEvery = 50,
      smoothWindow = NULL, log2Norm = FALSE, scaleTo = NULL)
    expect_identical(dim(result), c(1L, 2L))
    expect_equal(as.numeric(SummarizedExperiment::assay(result)), c(2.5, 6.5))
    expect_identical(rownames(result), "gene1")
  }
})

test_that("a single trajectory bin stays a feature by one matrix", {
  f <- trajectory_fixture(nfeatures = 3L)
  result <- GetTrajectory(f$object, assay = "RNA", groupEvery = 100,
    smoothWindow = NULL, log2Norm = FALSE, scaleTo = NULL)
  expect_identical(dim(result), c(3L, 1L))
  expect_equal(as.numeric(SummarizedExperiment::assay(result)), as.numeric(Matrix::rowMeans(f$counts)))
  expect_identical(rownames(result), rownames(f$counts))
  expect_identical(colnames(result), "T.0_100")
})

test_that("empty trajectory bins fail explicitly instead of propagating undefined means", {
  f <- trajectory_fixture(times = c(10,20,30,40))
  expect_error(GetTrajectory(f$object, assay = "RNA", groupEvery = 25), "empty trajectory bins")
  expect_error(GetTrajectory(f$object, assay = "RNA", groupEvery = 25,
    scaleTo = NULL, log2Norm = FALSE, smoothWindow = NULL), "empty trajectory bins")
  f$object$Trajectory <- rep(NA_real_, 4)
  expect_error(GetTrajectory(f$object, assay = "RNA", groupEvery = 25), "empty trajectory bins")
})

test_that("invalid trajectory metadata and absent TFs are diagnosed", {
  f <- trajectory_fixture()
  expect_error(GetTrajectory(f$object, assay = "RNA", trajectory.name = "absent"), "Cannot find")
  f$object$Trajectory <- rep("unknown", 8)
  expect_error(GetTrajectory(f$object, assay = "RNA"), "numeric")
  f$object$Trajectory <- seq(10, 80, 10)
  f$object$Trajectory[1] <- 101
  expect_error(GetTrajectory(f$object, assay = "RNA"), "between 0 and 100")
  a <- SummarizedExperiment::SummarizedExperiment(assays = list(mat = matrix(1:6, 1, dimnames = list("TF1", NULL))))
  b <- SummarizedExperiment::SummarizedExperiment(assays = list(mat = matrix(1:6, 1, dimnames = list("TF2", NULL))))
  expect_error(GetCorrelation(a,b), "Cannot find any shared features")
})

test_that("one shared TF yields the same correlation and p value as cor.test", {
  x <- c(1,3,2,5,4,8)
  y <- c(2,5,4,8,7,10)
  a <- SummarizedExperiment::SummarizedExperiment(assays = list(mat = rbind(TF1=x, other=6:1)))
  b <- SummarizedExperiment::SummarizedExperiment(assays = list(mat = rbind(unrelated=1:6, TF1=y)))
  result <- GetCorrelation(a,b)
  expected <- stats::cor.test(x,y)
  expect_identical(result$tfs, "TF1")
  expect_equal(result$correlation, unname(expected$estimate))
  expect_equal(result$p_value, expected$p.value)
  expect_equal(result$adj_p, expected$p.value)
})
