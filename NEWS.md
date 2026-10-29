## **scMEGA 1.1.0**

* Preserve RNA/ATAC count-to-cell correspondence in `CreatePairedObject()` and
  retain supported many-to-one mappings.
* Preserve singleton dimensions in trajectory and TF-correlation calculations.
* Validate Seurat v5 reference layers in `CoembedData()` while retaining legacy
  Assay and Signac ChromatinAssay support.
* Fix `TopEmbGRN()` to use its supplied graph and current igraph APIs.
* Declare runtime dependencies, add hosted package/security checks, and expand
  deterministic regression coverage.
* Update all three tutorials for current Seurat/Signac APIs. The Cardiomyocyte
  example uses a fixed, stratified 40,000-RNA/8,000-ATAC sample so exact CCA can
  run on a 64 GB workstation; it is not numerically identical to the full-data
  analysis.

Known limitations: ArchR remains a required GitHub dependency and four private
ArchR helpers are still used. Custom genomes and species-aware motif identifier
mapping remain separate feature/policy work.

## **scMEGA 0.2.0**

## **scMEGA 1.0.1**
