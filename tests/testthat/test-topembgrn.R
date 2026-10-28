test_that("TopEmbGRN uses the supplied graph for vertex types", {
  edges <- data.frame(
    from = c("A", "A", "A", "B", "B", "C", "E", "F", "B"),
    to = c("B", "C", "D", "C", "D", "D", "A", "E", "E")
  )
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE)
  igraph::V(graph)$type <- ifelse(
    igraph::V(graph)$name %in% c("A", "B"),
    "TF/Gene",
    "Gene"
  )

  expect_false(exists("netobj", envir = .GlobalEnv, inherits = FALSE))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  result <- TopEmbGRN(graph)

  expect_s3_class(result, "prcomp")
  expect_setequal(rownames(result$x), igraph::V(graph)$name)
  expect_equal(nrow(result$x), igraph::vcount(graph))
  expect_true(all(is.finite(result$x)))
})
