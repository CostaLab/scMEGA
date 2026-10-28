.ResolveAssayLayer <- function(object, assay, layer, caller = "This function") {
  if (is.null(assay) || !assay %in% SeuratObject::Assays(object)) {
    stop(caller, " requires an available assay; cannot find '", assay, "'.")
  }
  matches <- SeuratObject::Layers(object[[assay]], search = layer)
  if (length(matches) == 0L) {
    stop(caller, " requires layer '", layer, "' in assay '", assay, "'.")
  }
  if (length(matches) > 1L) {
    stop(
      caller, " requires one unambiguous '", layer, "' layer in assay '",
      assay, "', but found: ", paste(matches, collapse = ", "),
      ". Join only the intended disjoint layers with JoinLayers() before calling it."
    )
  }
  matches[[1L]]
}

.LayerDataExact <- function(object, assay, layer, cells = NULL,
                            features = NULL, caller = "This function") {
  resolved <- .ResolveAssayLayer(object, assay, layer, caller)
  SeuratObject::LayerData(
    object = object,
    assay = assay,
    layer = resolved,
    cells = cells,
    features = features
  )
}

.MotifNames <- function(object, assay, caller = "This function") {
  if (is.null(assay) || !assay %in% SeuratObject::Assays(object)) {
    stop(caller, " requires an available motif assay; cannot find '", assay, "'.")
  }
  assay.object <- object[[assay]]
  if (!inherits(assay.object, "ChromatinAssay")) {
    stop(caller, " requires assay '", assay, "' to be a ChromatinAssay with motifs.")
  }
  names <- unname(unlist(
    Signac::GetMotifData(assay.object, slot = "motif.names"),
    use.names = FALSE
  ))
  if (length(names) == 0L) {
    stop(caller, " requires motif names in assay '", assay, "'.")
  }
  names
}
