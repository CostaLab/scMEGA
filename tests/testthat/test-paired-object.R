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
  chromatin <- SeuratObject::CreateSeuratObject(
    Signac::CreateChromatinAssay(counts = f$atac[, rev(colnames(f$atac))],
      sep = c("-", "-"), min.cells = 10), assay = "ATAC")
  for (atac.object in list(f$obj.atac, chromatin)) {
    for (rows in list(seq_len(12), c(12L, seq_len(9)))) {
      pairs <- f$pairs[rows, ]
      result <- fn(pairs, obj.coembed = coembed, obj.rna = f$obj.rna,
        obj.atac = atac.object, rna.assay = "RNA", atac.assay = "ATAC")
      for (modality in c("rna", "atac")) {
        expected <- f[[modality]][, pairs[[toupper(modality)]], drop = FALSE]
        colnames(expected) <- pairs$cell_name
        expect_equal(as.matrix(SeuratObject::LayerData(result, assay = toupper(modality),
          layer = "counts")), as.matrix(expected))
      }
      expect_s4_class(result[["ATAC"]], "ChromatinAssay")
      expect_identical(result$source_id, setNames(pairs$RNA, pairs$cell_name))
      expected <- embedding[pairs$RNA, , drop = FALSE]
      rownames(expected) <- pairs$cell_name
      expect_equal(SeuratObject::Embeddings(result, "pca"), expected)
      expect_identical(as.integer(GenomicRanges::start(GenomicRanges::granges(result[["ATAC"]]))), c(100L, 300L))
    }
  }
})
