test_that("full ChromatinAssay output preserves counts, metadata, ranges and embeddings", {
  skip_if_not_installed("Signac")
  f <- pairing_fixture()
  coembed <- f$obj.rna
  coembed$source_id <- colnames(coembed)
  embedding <- matrix(seq_len(24), ncol = 2,
    dimnames = list(colnames(coembed), c("PC_1", "PC_2")))
  coembed[["pca"]] <- SeuratObject::CreateDimReducObject(embedding, key = "PC_", assay = "RNA")
  fn <- pairing_function()
  # Existing package lacks a Signac import: supply the actual constructor explicitly
  # in a child environment, without claiming this tests namespace installation.
  env <- new.env(parent = environment(fn))
  env$CreateChromatinAssay <- Signac::CreateChromatinAssay
  environment(fn) <- env
  result <- fn(f$pairs, obj.coembed = coembed, obj.rna = f$obj.rna,
    obj.atac = f$obj.atac, rna.assay = "RNA", atac.assay = "ATAC")
  for (modality in c("rna", "atac")) {
    expected <- f[[modality]][, f$pairs[[toupper(modality)]], drop = FALSE]
    colnames(expected) <- f$pairs$cell_name
    expect_equal(as.matrix(SeuratObject::LayerData(result, assay = toupper(modality),
      layer = "counts")), as.matrix(expected))
  }
  expect_s4_class(result[["ATAC"]], "ChromatinAssay")
  expect_identical(result$source_id, setNames(f$pairs$RNA, f$pairs$cell_name))
  expected <- embedding[f$pairs$RNA, , drop = FALSE]
  rownames(expected) <- f$pairs$cell_name
  expect_equal(SeuratObject::Embeddings(result, "pca"), expected)
  expect_identical(as.integer(GenomicRanges::start(GenomicRanges::granges(result[["ATAC"]]))), c(100L, 300L))
})
