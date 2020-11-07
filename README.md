# MitoX
MitoX: exploring mitochondrial heteroplasmies and gene expressions from single cell sequencing assays

## Introduction

MitoX is a computational framework to investigate mitochondrial heteroplasmies and mitochondrial gene expression in next generation sequencing data， including （sc)RNA-seq, (sc)ATAC-seq, across diverse platforms at cellular level. MitoX works on readily aligned (sc)RNA or (sc)ATAC-seq data. The MitoX functions gen_mut_profile() and gen_expr_profile() are the two main functions that takes as input either a single BAM/SAM file or a list of BAM/SAM files to generate the profiles for variants and gene expression. AnnData were used as the fundamental data structure which provides a scalable way of keeping track of data and learned annotations.


## Installation

**MitoX** depends on the following Python packages. These need to be installed separately:
```
pysam
HTSeq
anndata
```

To install **MitoX**, follow these instructions:

1. MitoX can be installed directly from Github using pip:

```alias
pip install git+https://github.com/huruifeng/MitoX.git
```

2. User can download the package from Github and install it locally:

```alias
git clone https://github.com/huruifeng/MitoX
cd MitoX
pip install .
```
