layer_fixture <- function(assay = c("assay5", "legacy", "chromatin"), split = FALSE) {
  assay <- match.arg(assay)
  cells <- paste0("cell", seq_len(20))
  counts <- Matrix::Matrix(
    matrix(seq_len(40), nrow = 2,
           dimnames = list(c("feature1", "feature2"), cells)),
    sparse = TRUE
  )
  assay.object <- switch(
    assay,
    assay5 = SeuratObject::CreateAssay5Object(
      counts = if (split) list(batch1 = counts[, 1:10], batch2 = counts[, 11:20]) else counts
    ),
    legacy = SeuratObject::CreateAssayObject(counts = counts),
    chromatin = Signac::CreateChromatinAssay(
      counts = `rownames<-`(counts, c("chr1-100-200", "chr2-300-400")),
      min.cells = 0,
      min.features = 0
    )
  )
  object <- SeuratObject::CreateSeuratObject(assay.object, assay = "RNA")
  object$Trajectory <- seq(5, 100, length.out = ncol(object))
  list(object = object, counts = counts)
}

test_that("trajectory layer access supports Seurat v5 assay classes without reordering", {
  results <- list()
  for (kind in c("assay5", "legacy", "chromatin")) {
    fixture <- layer_fixture(kind)
    results[[kind]] <- GetTrajectory(
      fixture$object,
      assay = "RNA",
      slot = "counts",
      groupEvery = 20,
      scaleTo = NULL,
      smoothWindow = NULL,
      log2Norm = FALSE
    )
    expect_identical(rownames(SummarizedExperiment::assay(results[[kind]])),
                     rownames(fixture$object))
  }
  expect_identical(SummarizedExperiment::assay(results$assay5),
                   SummarizedExperiment::assay(results$legacy))
  expect_identical(unname(SummarizedExperiment::assay(results$chromatin)),
                   unname(SummarizedExperiment::assay(results$legacy)))
})

test_that("ambiguous split layers fail until explicitly joined", {
  fixture <- layer_fixture("assay5", split = TRUE)
  expect_error(
    GetTrajectory(
      fixture$object,
      assay = "RNA",
      slot = "counts",
      groupEvery = 20,
      smoothWindow = NULL
    ),
    "one unambiguous 'counts' layer"
  )
  joined <- SeuratObject::JoinLayers(fixture$object, assay = "RNA")
  expect_no_error(
    GetTrajectory(
      joined,
      assay = "RNA",
      slot = "counts",
      groupEvery = 20,
      smoothWindow = NULL
    )
  )
})

test_that("paired-object construction rejects ambiguous split count layers", {
  fixture <- pairing_fixture()
  fixture$obj.rna[["RNA"]] <- SeuratObject::CreateAssay5Object(
    counts = list(
      batch1 = fixture$rna[, 1:6, drop = FALSE],
      batch2 = fixture$rna[, 7:12, drop = FALSE]
    )
  )
  expect_error(extract_pair_counts(fixture), "one unambiguous 'counts' layer")
  fixture$obj.rna <- SeuratObject::JoinLayers(fixture$obj.rna, assay = "RNA")
  actual <- extract_pair_counts(fixture)
  expect_identical(rownames(actual$rna), rownames(fixture$rna))
  expect_equal(
    as.matrix(actual$rna),
    as.matrix(`colnames<-`(
      fixture$rna[, fixture$pairs$RNA, drop = FALSE],
      fixture$pairs$cell_name
    ))
  )
})

test_that("public Signac motif accessor preserves motif identifiers", {
  fixture <- layer_fixture("chromatin")
  assay <- fixture$object[["RNA"]]
  motif.data <- Matrix::Matrix(
    matrix(c(1, 0, 0, 1), nrow = 2,
           dimnames = list(rownames(assay), c("motif1", "motif2"))),
    sparse = TRUE
  )
  motif.data <- methods::as(
    methods::as(motif.data, "generalMatrix"),
    "CsparseMatrix"
  )
  Signac::Motifs(assay) <- Signac::CreateMotifObject(
    data = motif.data,
    pwm = list(motif1 = matrix(1, 4, 1), motif2 = matrix(1, 4, 1)),
    motif.names = c(motif1 = "TF1", motif2 = "TF2")
  )
  fixture$object[["RNA"]] <- assay
  expect_identical(
    scMEGA:::.MotifNames(fixture$object, "RNA", "test"),
    c("TF1", "TF2")
  )
})
