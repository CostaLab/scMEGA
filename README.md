### scMEGA: Single-cell Multiomic Enhancer-based Gene regulAtory network inference

![](https://costalab.github.io/scMEGA/reference/figures/schematic.png)

You can install scMEGA via below commands:
```R
# Install devtools
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
# First install Seurate v5
devtools::install_github("satijalab/seurat", "seurat5", quiet = TRUE)

# Install Signac
devtools::install_github("stuart-lab/signac", "seurat5")
    
# Install scMEGA
devtools::install_github("CostaLab/scMEGA")
```

We provided the following tutorials to show how to use scMEGA to build GRN by using single-cell multiomics/multimodal data:

* [Myofibroblasts differentiation](https://costalab.github.io/scMEGA/articles/myofibroblast-GRN.html)

* [Cardiomyocytes remodelling](https://costalab.github.io/scMEGA/articles/cardiomyocyte-GRN.html)

* [CD4 T cells activcation using 10X multiome data](https://costalab.github.io/scMEGA/articles/pbmc_10x_multiome.html)

Please consider citing our paper if you used scMEGA:
```
@article{li2023scmega,
  title={scMEGA: single-cell multi-omic enhancer-based gene regulatory network inference},
  author={Li, Zhijian and Nagai, James S and Kuppe, Christoph and Kramann, Rafael and Costa, Ivan G},
  journal={Bioinformatics Advances},
  volume={3},
  number={1},
  pages={vbad003},
  year={2023},
  publisher={Oxford University Press}
}
```

