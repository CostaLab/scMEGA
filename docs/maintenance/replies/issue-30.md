Thank you for the follow-up. A GRN inferred along the naïve-to-memory CD4 trajectory describes associations in the cells/states included in that analysis. It should not be described as a GRN for all PBMC cell types.

In the current implementation, `SelectTFs()`, `SelectGenes()` and `GetTFGeneCorrelation()` all call `GetTrajectory()`. `trajectory.name` selects an existing numeric metadata column; it is not a switch to disable trajectory analysis. Renaming it to `"All"`, setting it to NULL, or assigning arbitrary cell-order values does not provide a valid trajectory-free analysis.

You can inspect participation without modifying the object:

```r
t <- pbmc.all.cells[["Trajectory"]][, 1]
table(is.na(t))
sum(t > 0 & t <= 100, na.rm=TRUE)  # cells included by current default interval bins
```

Having all PBMC cells in the object does not mean they all contributed: missing trajectory values are excluded. The current method aggregates and smooths values over trajectory bins before computing associations. For heterogeneous PBMC populations, choose biologically justified lineages/states and report their networks separately, or design a separately validated trajectory-free method. There is no existing parameter change that establishes that new statistical method.

The inferred edges combine motif evidence, selected peak-to-gene associations and TF-activity/gene-expression correlations. They are computational predictions, not experimentally established causal regulation. Correlation/FDR values describe the processed-bin association calculation; `n_peaks` counts supplied motif/peak-link support rows, not a probability of regulation. `GetGRN()` assembles supplied tables and does not itself make the analysis trajectory-free or automatically remove missing/negative correlations; input selection/filtering matters.

Verification: the function call chain was inspected, and a synthetic `GetGRN()` example verified its joins, support counts and retention of missing/negative correlations. No all-PBMC biological analysis was performed and no source change was made. Supporting trajectory-free inference would require a separate methodological proposal and validation, rather than bypassing the trajectory checks.
