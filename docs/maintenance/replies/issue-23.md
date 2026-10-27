Thank you for the traceback and follow-up details. The shown failure happens while retrieving an assay layer, before genome-specific peak-to-gene annotation. It does not by itself identify a rat-genome problem.

Your traceback uses an older `GetTrajectory()` implementation with `GetAssayData()`. Current checked source uses `LayerData()`, but the requested data must still exist. `SelectTFs()` requests `data` layers from the RNA and chromVAR assays. Converting an assay class does not create the missing normalized expression or chromVAR activity.

Please inspect the actual object and installed code:

```r
body(scMEGA::GetTrajectory)
for (a in c("SoupXRNA", "chromvar", "ATAC")) {
  print(a)
  print(class(objG[[a]]))
  print(SeuratObject::Layers(objG[[a]], search=NA))
}
Signac::GetMotifData(objG, assay="ATAC", slot="motif.names")
sessionInfo()
```

Use your actual RNA assay name. If `SoupXRNA` contains the intended corrected counts but no normalized data, complete the appropriate RNA normalization, for example `Seurat::NormalizeData(objG, assay="SoupXRNA")`, retaining its returned object. If chromVAR activity is absent, complete motif annotation and `RunChromVAR()` on the correct object. Do not replace TF activity with RNA counts or normalize chromVAR deviation scores as if they were RNA counts.

Keep the peak assay as `ChromatinAssay` so its genomic ranges/motifs remain available. The chromVAR activity assay itself can be an ordinary `Assay` with a `data` layer; it does not need to be cast to ChromatinAssay. Once inputs are valid, use the existing argument instead of copying/modifying the function to change the RNA name:

```r
sel <- SelectTFs(
  object=objG, rna.assay="SoupXRNA", tf.assay="chromvar", atac.assay="ATAC",
  trajectory.name="Trajectory", cor.cutoff=0.1
)
```

Verification: a synthetic counts-only RNA assay reproduced a missing-data failure, and appropriate RNA normalization restored `GetTrajectory()` output without casting. A data-only ordinary assay also produced correct TF-activity bin means. Full `SelectTFs()` on your rat data was not tested. Motif-to-gene naming (#18/#33), singleton dimensions (#25), split layers (#27) and package imports are separate issues. No blanket downgrade or historical development-branch installation is recommended; please share the inspection output and a small synthetic or anonymized reproduction if the error persists. Local fixes/proposals are not yet published.
