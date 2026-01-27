FROM rocker/r-ver:4.4.1

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    pkg-config \
    zlib1g-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
  && rm -rf /var/lib/apt/lists/*

RUN R -e 'options(repos=c(CRAN="https://cloud.r-project.org")); \
          install.packages(c("data.table","optparse","ggplot2","BiocManager")); \
          options(repos=BiocManager::repositories()); \
          BiocManager::install("DESeq2", ask=FALSE, update=FALSE); \
          stopifnot(requireNamespace("DESeq2", quietly=TRUE)); \
          stopifnot(requireNamespace("data.table", quietly=TRUE))'