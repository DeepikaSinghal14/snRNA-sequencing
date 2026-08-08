# R environment for the snRNA-seq pipeline.
# Build: docker build -t ghcr.io/deepikasinghal14/snrna-seq-r:0.1.0 -f environment/Dockerfile .
FROM bioconductor/bioconductor_docker:RELEASE_3_18

RUN R -e "install.packages(c('optparse','dplyr','tidyverse','data.table','openxlsx','readxl','rio','patchwork','cowplot','viridis','qs'), repos = 'https://cloud.r-project.org')" \
    && R -e "BiocManager::install(c('SingleCellExperiment','DESeq2','scater','biomaRt'), update = FALSE, ask = FALSE)" \
    && R -e "install.packages('Seurat', repos = 'https://cloud.r-project.org')" \
    && R -e "install.packages('remotes', repos = 'https://cloud.r-project.org'); remotes::install_github(c('constantAmateur/SoupX', 'samuel-marsh/scCustomize'))"

WORKDIR /workspace
