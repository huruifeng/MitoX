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

```
pip install git+https://github.com/huruifeng/MitoX.git
```

2. User can download the package from Github and install it locally:

```
git clone https://github.com/huruifeng/MitoX
cd MitoX
pip install .
```
## Main functions
The full documents for all the functions in MitoX can be accessed in the [**Wiki page**](https://github.com/huruifeng/MitoX/wiki/MitoX)

1. Generate the **heteroplamy profiles** for all samples/cells
```
 mut_adata = gen_mut_profile(bams, fmt="bam", sp="human", chr_name="MT", seq_type="sc", combined=False, tag_name="CB",
                           barcodes="NULL", min_read=200, max_depth=1e5, min_baseq=25, min_mapq=0, n_jobs=2)
    
 * Parameters:
    - bams      (str): path(s) pointing to BAM alignment file(s).
    - fmt       (str): File type, bam or sam.
    - sp        (str): Species, MitoX supports 'human' and 'mus'.
    - chr_name  (str): Name of mitochondrial genome as specified in the BAM files.
    - seq_type  (str): Sequencing type:bulk or sc. MitoX supports both bulk sequencing data and single cell sequencing data.
    - combined  (str): If the BAM is merged file or not (It should be set to True for 10X or droplet data).
    - tag_name  (str): The name of the tag corresponding the cellular barcode. Default = "CB". For droplet scRNA-seq only.
    - barcodes (list): The barcode list corresponding to the cells. For 10X genomics scRNA-seq data only.
    - min_read  (int): The minimum number of read counts to be considered a valid barcode (cell) in the analysis. Default = 200. 
                       For droplet scRNAseq technologies only.
    - max_depth (int): The maximum depth of reads considered at any position.
    - min_baseq (int): The minimum read base quality below which the base is ignored.
    - min_mapq  (int): The minimum map quality below which the read is ignored.
    - n_jobs    (int): The number of threads will be used,each thread can handle one BAM file.

*return (AnnData): anndata object containd the heteroplamy profiles
```
2. Generate the mitochondrial **gene expression profiles** for all samples/cells
```
 expr_adata = gen_expression_profile(bams, gtf_file, fmt="bam", sp="human", chr_name="MT", eq_type="sc",
                                    combined=False, tag_name="CB", barcodes="NULL", min_read=200, stranded=True,
                                    feature_type="exon", min_mapq=10, n_jobs=2):
    
 * Parameters:
    - bams      (str): path(s) pointing to BAM alignment file(s).
    - gtf_file  (str): Location of the GTF file that will be used for gene quantification. 
    - fmt       (str): File type, bam or sam.
    - sp        (str): Species, MitoX supports 'human' and 'mus'.
    - chr_name  (str): Name of mitochondrial genome as specified in the BAM files.
    - seq_type  (str): Sequencing type:bulk or sc. MitoX supports both bulk sequencing data and single cell sequencing data.
    - combined  (str): If the BAM is merged file or not (It should be set to True for 10X or droplet data).
    - tag_name  (str): The name of the tag corresponding the cellular barcode. Default = "CB". For droplet scRNA-seq only.
    - barcodes (list): The barcode list corresponding to the cells. For 10X genomics scRNA-seq data only.
    - min_read  (int): The minimum number of read counts to be considered a valid barcode (cell) in the analysis. Default = 200. 
                       For droplet scRNAseq technologies only.
    - stranded (bool): Whether the data is from a strand-specific assay (default: yes).
                       For stranded=no, a read is considered overlapping with a feature regardless of whether it is mapped 
                       to the same or the opposite strand as the feature. 
                       For stranded=yes and single-end reads, the read has to be mapped to the same strand as the feature. 
                       For paired-end reads, the first read has to be on the same strand and the second read on the opposite strand. 
                       
    - min_mapq  (int): The minimum map quality below which the read is ignored.
    - n_jobs    (int): The number of threads will be used,each thread can handel one BAM file.
    - feature_type (str): Feature type (3rd column in GTF file) to be used. 
                          All features of other type are ignored. (default, suitable for RNA-Seq analysis using an Ensembl GTF file: exon).

*return (AnnData): anndata object containd the gene expression profiles
```

## Examples

* Example for scRNA-seq: [1.1-STREAM_scRNA-seq (Bifurcation).ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/1.1.STREAM_scRNA-seq%20%28Bifurcation%29.ipynb?flush_cache=true)

* Example for scRNA-seq: [1.2-STREAM_scRNA-seq (Multifurcation).ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/1.2.STREAM_scRNA-seq%20%28Multifurcation%29.ipynb?flush_cache=true)

* Example for scATAC-seq(using peaks): [2.1-STREAM_scATAC-seq_peaks.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/2.1-STREAM_scATAC-seq_peaks.ipynb?flush_cache=true)

* Example for scATAC-seq(using k-mers): [2.2-STREAM_scATAC-seq_k-mers.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/2.2.STREAM_scATAC-seq_k-mers.ipynb?flush_cache=true)

* Example for scATAC-seq(using motifs): [2.3-STREAM_scATAC-seq_motifs.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/2.3.STREAM_scATAC-seq_motifs.ipynb?flush_cache=true)

* Example for *mapping* feature: [3-STREAM_mapping.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/3.STREAM_mapping.ipynb?flush_cache=true)

* Example for complex trajectories: [4-STREAM_complex_trajectories.ipynb](https://nbviewer.jupyter.org/github/pinellolab/STREAM/blob/master/tutorial/4.STREAM_complex_trajectories.ipynb?flush_cache=true)

## Document and API
The full descriptive documents for all the functions and APIs in MitoX can be accessed in the [**Wiki page**](https://github.com/huruifeng/MitoX/wiki/MitoX)

## Citaion
## Note

