#!/bin/bash

#R -e "devtools::document()"
#R -e "pkgdown::build_site(preview = FALSE, lazy=TRUE)"
R -e "pkgdown::build_article(name='pbmc_10x_multiome', lazy=TRUE)"
