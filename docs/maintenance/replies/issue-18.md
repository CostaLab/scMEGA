Thank you for reporting the mouse-symbol mismatch and sharing the workaround. The source confirms the limitation: motif labels are assigned to activity rows and then matched exactly to RNA feature names. A case change can recover simple examples, but it is not a general species-aware identifier mapping.

Please retain motif IDs and a mapping table rather than replacing labels globally. For example:

```r
motif_names <- unlist(Signac::GetMotifData(object, assay="ATAC", slot="motif.names"))
activity <- SeuratObject::LayerData(object, assay="chromvar", layer="data")
rna <- SeuratObject::LayerData(object, assay="RNA", layer="data")
ids <- rownames(activity)
label <- unname(motif_names[ids])
data.frame(motif_id=ids, motif_label=label, exact_RNA_match=label %in% rownames(rna))
```

Changing case may collapse distinct labels. Composite motifs are not equivalent to independent monomer TFs, and multiple motif IDs can represent the same TF or family. Signac can also append suffixes to make duplicate names unique. Those suffixes are not gene symbols. A curated mapping must identify the species and annotation release and define how aliases, complexes and multiple motifs per gene are handled.

Apply any reviewed mapping consistently to TF selection, `GetTFGeneCorrelation()`, and the TF identifiers used in `GetGRN()`'s motif-matching matrix. Changing only one function can cause missing-index errors or mismatched network edges. Rerun affected selection/correlation steps: changing the set of matched TFs also changes the multiple-testing universe. Do not silently discard missing TFs or duplicate motif evidence.

Verification: synthetic inputs with real Signac motif metadata reproduced exact-name losses and case/composite ambiguities. This supports a dedicated mapping enhancement, not an unconditional title-case patch. No species-mapping code was changed, and no PR has been published. Please provide representative motif-ID/name/RNA-ID examples and database provenance so the mapping policy can be reviewed. The later Seurat missing-layer discussion is a separate compatibility problem addressed in the draft response for #23; it is not solved by renaming genes.
