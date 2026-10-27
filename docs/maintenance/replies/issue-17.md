Thank you for the question and the custom-annotation examples. The unmodified package currently supports only `hg19`, `hg38`, `mm9` and `mm10`. `SelectGenes()` passes its `genome` argument to `PeakToGene()`, which selects one of the four bundled annotations. It does not currently accept a custom GRanges/TxDb annotation through that argument.

Adding another species is feasible, but needs a reviewed annotation interface and validation of assembly, chromosome names, coordinate convention, strand/TSS definitions, gene IDs/symbols, duplicate transcripts and matrix-row alignment. Renaming `chr1` to `1` can reconcile naming for the **same assembly**; it cannot convert coordinates between assemblies or species. Do not substitute the nearest supported genome.

The posted macaque `SelectGenes()` replacement selects an annotation object but then constructs placeholder `tfs` values and a random heatmap; it does not implement peak-to-gene inference. The annotation selection would need to reach a generalized `PeakToGene()` implementation, preserving its actual statistical calculations and filters. The user-modified rat example is useful input for feature design, but it does not establish general validated support.

The later `tf_activity[tf.use, ]` error is a separate problem: at least one requested TF name is absent from the processed activity matrix (or its dimensions have changed). Check consistent motif naming and the same object/trajectory parameters across TF selection and correlation:

```r
requested <- df.cor$tfs
motif_names <- unname(unlist(Signac::GetMotifData(
  pbmc.t.cells, assay="ATAC", slot="motif.names"
)))
setdiff(requested, motif_names)
```

This is a diagnostic, not permission to silently remove the missing TFs. If that difference is empty, please provide `tf.use`, `rownames(tf_activity)` at the failing line, the actual function body being used, and `sessionInfo()`; internal filtering can also remove activity rows.

Verification: a synthetic call to `PeakToGene(genome="rn6")` reproduces the unsupported-genome error. Three synthetic links using the bundled mm10 annotation and actual C++ correlation routine preserve peak/gene identity and match direct correlations after gene-row reordering. No rat/macaque dataset or custom-genome extension was validated. Please provide the desired species/assembly accession, annotation source/release, RNA identifier scheme, and example peak coordinates before a separate feature PR is designed. Also specify the genome explicitly: `SelectGenes()` defaults to hg38, while direct `PeakToGene()` defaults to hg19.
