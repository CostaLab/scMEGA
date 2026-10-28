# Load only the relevant source for focused tests when package installation is blocked.
# Normal R CMD check uses the installed package function.
pairing_function <- function() {
  if (Sys.getenv("SCMEGA_SOURCE_TESTS") == "true") {
    env <- new.env(parent = asNamespace("SeuratObject"))
    sys.source(file.path(Sys.getenv("SCMEGA_SOURCE_ROOT"), "R", "pair_cells.R"), env)
    env$CreatePairedObject
  } else {
    getExportedValue("scMEGA", "CreatePairedObject")
  }
}

pairing_fixture <- function(assay5 = TRUE, sparse = TRUE) {
  cells <- paste0("c", seq_len(12))
  rna <- matrix(seq_len(24), 2, dimnames = list(c("geneB", "geneA"), cells))
  atac <- matrix(101:124, 2,
    dimnames = list(c("chr1-100-200", "chr2-300-400"), cells))
  if (sparse) {
    rna <- Matrix::Matrix(rna, sparse = TRUE)
    atac <- Matrix::Matrix(atac, sparse = TRUE)
  }
  make <- if (assay5) SeuratObject::CreateAssay5Object else SeuratObject::CreateAssayObject
  list(rna = rna, atac = atac,
    obj.rna = SeuratObject::CreateSeuratObject(make(counts = rna), assay = "RNA"),
    obj.atac = SeuratObject::CreateSeuratObject(make(counts = atac[, rev(cells)]), assay = "ATAC"),
    pairs = data.frame(RNA = rev(cells), ATAC = rev(cells),
      cell_name = paste0("pair", seq_along(cells))))
}

# Stop immediately before object construction, after the real extraction/rename code.
# Only the constructor is replaced; LayerData dispatch and sparse operations are real.
extract_pair_counts <- function(fixture, pairs = fixture$pairs) {
  fn <- pairing_function()
  env <- new.env(parent = environment(fn))
  env$CreateSeuratObject <- function(counts, ...) {
    stop(structure(list(message = "capture", call = NULL, rna = counts,
      atac = get("atac.counts", envir = parent.frame())),
      class = c("pair_counts", "error", "condition")))
  }
  environment(fn) <- env
  tryCatch(fn(pairs, obj.rna = fixture$obj.rna, obj.atac = fixture$obj.atac,
    rna.assay = "RNA", atac.assay = "ATAC"), pair_counts = identity)
}
