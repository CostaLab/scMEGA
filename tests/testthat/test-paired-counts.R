test_that("pair identity and feature counts survive assay ordering and representation", {
  for (assay5 in c(TRUE, FALSE)) {
    for (sparse in c(TRUE, FALSE)) {
      f <- pairing_fixture(assay5, sparse)
      for (order in list(seq_len(12), rev(seq_len(12)), c(3L, 1L, 7L))) {
        pairs <- f$pairs[order, ]
        actual <- extract_pair_counts(f, pairs)
        for (modality in c("rna", "atac")) {
          expected <- f[[modality]][, pairs[[toupper(modality)]], drop = FALSE]
          colnames(expected) <- pairs$cell_name
          expect_equal(as.matrix(actual[[modality]]), as.matrix(expected))
          expect_identical(rownames(actual[[modality]]), rownames(expected))
          expect_s4_class(actual[[modality]], "dgCMatrix")
        }
      }
    }
  }
})

test_that("declared pair mappings allow distinct barcodes and repeated source cells", {
  f <- pairing_fixture()
  pairs <- f$pairs[c(1, 1, 2), ]
  pairs$ATAC <- c("c2", "c3", "c4")
  pairs$cell_name <- c("p1", "p2", "p3")
  actual <- extract_pair_counts(f, pairs)
  expected <- f$rna[, pairs$RNA, drop = FALSE]
  colnames(expected) <- pairs$cell_name
  expect_equal(as.matrix(actual$rna), as.matrix(expected))
  expected <- f$atac[, pairs$ATAC, drop = FALSE]
  colnames(expected) <- pairs$cell_name
  expect_equal(as.matrix(actual$atac), as.matrix(expected))
  single <- extract_pair_counts(f, pairs[1, , drop = FALSE])
  expect_identical(dim(single$rna), c(2L, 1L))
  expect_identical(dim(single$atac), c(2L, 1L))
})

test_that("invalid pair identifiers fail explicitly instead of silently reassigning counts", {
  f <- pairing_fixture()
  for (column in c("RNA", "ATAC", "cell_name")) {
    for (bad in list(NA_character_, "")) {
      pairs <- f$pairs
      pairs[[column]][1] <- bad
      expect_error(extract_pair_counts(f, pairs), "identifiers")
    }
  }
  pairs <- f$pairs
  pairs$cell_name[2] <- pairs$cell_name[1]
  expect_error(extract_pair_counts(f, pairs), "unique")
  expect_error(extract_pair_counts(f, f$pairs[FALSE, ]), "non-empty")
  expect_error(extract_pair_counts(f, f$pairs[, c("RNA", "ATAC")]), "cell_name")
  for (column in c("RNA", "ATAC")) {
    pairs <- f$pairs
    pairs[[column]][1] <- "missing-cell"
    expect_error(extract_pair_counts(f, pairs), "counts")
  }
})


test_that("factor IDs use names and numeric positions are rejected", {
  f <- pairing_fixture()
  pairs <- as.data.frame(lapply(f$pairs, factor))
  actual <- extract_pair_counts(f, pairs)
  expected <- f$rna[, as.character(pairs$RNA), drop = FALSE]
  colnames(expected) <- as.character(pairs$cell_name)
  expect_equal(as.matrix(actual$rna), as.matrix(expected))
  pairs <- f$pairs
  pairs$RNA <- seq_len(nrow(pairs))
  expect_error(extract_pair_counts(f, pairs), "identifiers")
})

test_that("aligned and partially overlapping source cell sets preserve declared pairs", {
  f <- pairing_fixture()
  for (atac.cells in list(colnames(f$atac), colnames(f$atac)[-1])) {
    f$obj.atac <- SeuratObject::CreateSeuratObject(
      SeuratObject::CreateAssay5Object(counts = f$atac[, atac.cells, drop = FALSE]),
      assay = "ATAC")
    pairs <- f$pairs[rev(seq_len(nrow(f$pairs))), ]
    pairs <- pairs[pairs$ATAC %in% atac.cells, ]
    actual <- extract_pair_counts(f, pairs)
    for (modality in c("rna", "atac")) {
      expected <- f[[modality]][, pairs[[toupper(modality)]], drop = FALSE]
      colnames(expected) <- pairs$cell_name
      expect_equal(as.matrix(actual[[modality]]), as.matrix(expected))
    }
  }
})
