# MitoX
MitoX: exploring mitochondrial heteroplasmies and gene expressions from single cell sequencing assays

## Introduction

MitoX is a python package for the analysis of mitochondrial variants from both bulk and single cell RNA-seq and ATAC-seq data. 
The two main function ***gen_mut_profile()*** and ***gen_expr_profile()*** generate the mitochondrial heteroplasmy profiles 
and mitochondrial gene expression profiles. After generating the heteroplasmy profiles or gene expression profiles, 
MitoX provides various analysis and powerful plot functions. 
One of the major functions is samples clustering, MitoX implements commonly used dimension reduction methods,
 such as `PCA`, `t-SNE`, and `UMAP` through the ***dimension_reduction()*** function. 
 User can do samples clustering based on the results from ***dimension_reduction()*** using ***plot_clusters()**** . 
 In plot_clusters(), MitoX can dtetcte samles clusters automatically using the ***Density-Based Spatial Clustering of Applications with Noise (DBSCAN)***
 or ***KMeans*** algorithm, and color the clusters. 
 Other major funtions include ***plot a clustered heatmap*** among variants and samples, 
 ***plot read coverage*** across the mitochondrial whole genome, ***filter low occurrence variants***, ***filter samples with few variants***, 
 ***calculate the genotype similarity*** based the heteroplasmy profiles, ***compare gene expression levels*** among user defined groups, 
 or ***plot gene expression profiles among clustered groups***. 

Our MitoX computational framework supports exploring mitochondrial variants and gene expressions from various sequencing assays for **both human and mouse mitochondrion**. 
MitoX can be applied to both of **well-based** and **droplet-based** single cell sequencing data formats and thus provides universal compatibility with all major single cell sequencing platforms. 

MitoX also support multiple threads to speed up the calculations. The typical run time on sequencing data with BAM file size of 1 GB is less than 30s on a 2.4 GHz CPU with less than 1 GB of memory. 
Utilizing the supports of multiple threads running, MitoX can process several BAM files at the same time, 
so **MitoX is a highly efficient tool for analyzing mitochondrial heteroplasmy and gene expressions.**

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
The full documents for all functions in MitoX can be accessed in the [**Wiki page**](https://github.com/huruifeng/MitoX/wiki/MitoX)

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

* Example 1: [Deatcte heteroplasmy profiles from (sc)RNA-seq and (sc)ATAC-seq data]()

* Example 2: [Investigate mitochondrial gene expressions and plot results]()

* Example 3: [Running MitoX on droplet based dataset (10X dataset)]()

## Document and API
1. The full descriptive documents for MitoX can be accessed on the [**Wiki page**](https://github.com/huruifeng/MitoX/wiki/MitoX)
2. User can also open the [**Doc folder**](https://github.com/huruifeng/MitoX/tree/master/Doc) for API reference.

## Citaion


