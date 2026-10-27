Thank you for the report. The screenshot shows an argument-name mismatch: `object`, `use.assay1` and `use.assay2` are not parameters of the current `CreatePairedObject()` implementation. The old example used those names.

The current call is:

```r
obj.pair <- CreatePairedObject(
  df.pair = df.pair,
  obj.coembed = obj.coembed2,
  obj.rna = obj.rna,
  obj.atac = obj.atac,
  rna.assay = "RNA",
  atac.assay = "ATAC"
)
```

Use the actual assay names in your separate RNA and ATAC objects (for example, `"peaks"` if that is your accessibility assay). `obj.coembed` supplies RNA-cell metadata and reductions; `obj.rna` and `obj.atac` supply the original counts. `df.pair` needs the columns `RNA`, `ATAC` and `cell_name`, with source cell identifiers and unique output names.

Verification: the old call reproduces the screenshot's `unused arguments` error. The current signature and real paired-object call were verified in the local regression suite. The directly related function documentation has been corrected locally. The separate #32 count-order correction also remains locally verified; no PR or release containing these changes has been published. The final isolated environment supports installation; with the separate metadata proposal, installed tests pass, but R CMD check still reports six warnings and three notes. This diagnoses the screenshot error; any later runtime failure would need the exact call, assay layers/classes and `sessionInfo()`.
