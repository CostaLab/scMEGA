Thank you for confirming that the component approach produced a plot. `NetCentPlot()` uses a focus layout that requires a connected graph. A general network layout can display disconnected components, so the two plots need not behave the same way.

For your follow-up, use the existing `highlights` argument to label the chosen node:

```r
netobj <- igraph::graph_from_data_frame(df.grn2, directed=TRUE)
igraph::V(netobj)$type <- ifelse(
  igraph::V(netobj)$name %in% df.grn2$tf, "TF/Gene", "Gene"
)
focus <- "RUNX1"
i <- match(focus, igraph::V(netobj)$name)
stopifnot(!is.na(i))
cc <- igraph::components(netobj, mode="weak")
net.focus <- igraph::induced_subgraph(
  netobj, which(cc$membership == cc$membership[i])
)
p <- NetCentPlot(net.focus, gene=focus, highlights=focus)
p
```

Construct the graph before assigning its vertex attributes. Match vertex **names**, rather than comparing numeric vertex positions with `"RUNX1"`. The default labeling branch only labels nodes with `type == "TF/Gene"` and betweenness above the graph mean; assigning the type `"TF"` will not satisfy that condition. `highlights="RUNX1"` selects its label explicitly.

Verification: a synthetic disconnected graph reproduced the focus-layout error. Selecting the weakly connected component and supplying `highlights` produced a rendered plot with RUNX1 labeled; the original graph was unchanged. Your attached plot was inspected. No package source change is needed. Component subsetting restricts the display; it does not recompute the inferred GRN or make the other components biologically irrelevant.
