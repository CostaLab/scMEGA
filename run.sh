#!/bin/bash

# R -e "devtools::document()"
# R -e "pkgdown::build_site(preview = FALSE, lazy=TRUE)"
# R -e "pkgdown::build_article(name='articles/myofibroblast-GRN', quiet=FALSE)"
R -e "pkgdown::build_article(name='articles/install', quiet=FALSE)"
R -e "pkgdown::build_article(name='articles/pbmc_10x_multiome', quiet=FALSE)"
