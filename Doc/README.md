.markdown-block {
    background: #101010;
}

# MitoX
MitoX: exploring mitochondrial heteroplasmies and gene expressions from single cell sequencing assays

Table of Content
================
* [Overview](#overview)
* [Note (important !!!)](#note)
* [Functions & API](#api)
    * [gen_mut_profile()](#gen_mut_profile)
    * [gen_expression_profile()](#gen_expression_profile)
    * [read_sam()](#read_sam)
    * [read_table()](#read_table)
    * [write_table()](#write_table)
    * [get_coverage()](#get_coverage)
    * [plot_coverage()](#plot_coverage)
    * [add_metadata()](#add_metadata)
    * [filter_features()](#filter_features)
    * [filter_cells()](#filter_cells)
    * [gen_expression_profile](#gen_expression_profile)
    * [gen_expression_profile](#gen_expression_profile)
    * [gen_expression_profile](#gen_expression_profile)
    * [gen_expression_profile](#gen_expression_profile)
    * [gen_expression_profile](#gen_expression_profile)


    

<a name="overview"/>

## Overview

MitoX is a computational framework to investigate mitochondrial heteroplasmies and mitochondrial gene expression in next generation sequencing data， including (sc)RNA-seq, (sc)ATAC-seq, across diverse platforms at cellular level. MitoX works on readily aligned (sc)RNA or (sc)ATAC-seq data. The MitoX functions gen_mut_profile() and gen_expr_profile() are the two main functions that takes as input either a single BAM/SAM file or a list of BAM/SAM files to generate the profiles for variants and gene expression. AnnData were used as the fundamental data structure which provides a scalable way of keeping track of data and learned annotations.

<a name="note"/>

## Note (important !!!)
* The BAM/SAM file name will be used as sample or cell names.
* In MitoX, when preparing dataframe, sample names (cells) should be in columns, and features (mutations,gene names) are in rows.
* When using add_metadata() function, the first column must be the sample names, 
and the second column MUST be the **'label'** column to indicate the sample groups.
* In function **combine_mut_expr(mut_adata, expr_adata)**, the first parameter must be mutation anndata object.


<a name="gen_mut_profile"/>

## Functions & API

> 1. Generate the **heteroplamy profiles** for all samples/cells
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

*return (AnnData): anndata object that contains the heteroplamy profiles
```
<a name="gen_expression_profile"/>

> 2. Generate the mitochondrial **gene expression profiles** for all samples/cells
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

*return (AnnData): anndata object that contains the gene expression profiles
```

<a name="read_sam"/>

> 3. Read SAM/BAM files and extract coverage data for all samples/cells. This function has been integrated in to **gen_mut_profile()**, 
>but user can still ues it, the returned results can be the input of gen_mut_profile().
```
 samples = read_sam(bam_folder, fmt="bam",sp="human", chr_name="MT", seq_type="sc",combined=False, tag_name="CB", barcodes="NULL",
                       min_read=200, min_baseq=25, min_mapq=0, n_jobs=2):
    
 * Parameters:
    - bam_folder(str): path(s) pointing to BAM alignment file(s).
    - fmt       (str): File type, bam or sam.
    - sp        (str): Species, MitoX supports 'human' and 'mus'.
    - chr_name  (str): Name of mitochondrial genome as specified in the BAM files.
    - seq_type  (str): Sequencing type:bulk or sc. MitoX supports both bulk sequencing data and single cell sequencing data.
    - combined  (str): If the BAM is merged file or not (It should be set to True for 10X or droplet data).
    - tag_name  (str): The name of the tag corresponding the cellular barcode. Default = "CB". For droplet scRNA-seq only.
    - barcodes (list): The barcode list corresponding to the cells. For 10X genomics scRNA-seq data only.
    - min_read  (int): The minimum number of read counts to be considered a valid barcode (cell) in the analysis. Default = 200. 
                       For droplet scRNAseq technologies only.
    - min_baseq (int):
    - min_mapq  (int): The minimum map quality below which the read is ignored.
    - n_jobs    (int): The number of threads will be used,each thread can handel one BAM file.

*return (Sample): a sample or sample list that contiains coverage information for each sample.
```

<a name="read_table"/>

> 4. Read data from text file.
```
 x_adata = read_table(file,sep="\t",type="mut"):
    
 * Parameters:
    - file  (str): path pointing to file.
    - sep   (str): Delimiter to use. 
    - type  (str): Data type - mutation profile table ('mut') or gene expression profile table ('expr').

*return (AnnData): anndata object that contains the mutation / gene expression profiles
```
<a name="write_table"/>

> 5. Write AnnData into local file.
```
write_table(adata,file,sep="\t"):
    
 * Parameters:
    - adata (AnnData): anndata object that contains the mutation / gene expression profiles.
    - file      (str): file path and name 
    - sep       (str): String of length 1. Field delimiter for the output file..

*return : NOne
```

<a name="get_coverage"/>

> 6. Get the coverage information for sample(s).
```
cov = get_coverage(adata,name="all")
    
 * Parameters:
    - adata (AnnData): anndata object that contains the mutation / gene expression profiles.
    - name      (str): sample name, when "all", returns a dataframe for all samples 

*return (Series or DataFrame): coverage data for sample(s) 
```
<a name="plot_coverage"/>

> 7. Plot the coverage information for sample(s).
```
plot_coverage(x_df,color=None,fig_name="coverage_plot.png",fig_size=None,plot_type="circle",log2=True)
    
 * Parameters:
    - x_df (Series,DataFrame): A Series of DataFrame that contains the coverage data. 
                               It can be the return value from get_coverage().
    - color       (list, str): A color name (e.g. 'red') for the Series，
                               OR, a list of colors for the DataFrome.
    - fig_name          (str): figure loaction and name for saving.
    - fig_size (float, float): Width, height in inches.
    - plot_type         (str): "circle" or "plane".
    - log2             (bool): Do log2 transform on coverage values or not.

*return:  None 
```
<a name="add_metadata"/>

> 8. Add annotation data (metadata) to samples.
```
add_metadata(adata,file_name,delimiter='\t')
    
 * Parameters:
    - adata  (AnnData): anndata object that contains the mutation / gene expression profiles.
    - file_name  (str): file path and name to the metadata file. 
    - sep        (str): String of length 1. Field delimiter of the metadata file.

*return:  None 
```
<a name="filter_features"/>

> 9. Filter out features based on different metrics.
```
filter_features(adata, min_n_cells = None, max_n_cells=None, min_pct_cells = None, max_pct_cells=None)
    
 * Parameters:
    - adata           (AnnData): anndata object that contains the mutation / gene expression profiles.
    - min/max_n_cells     (int): Minimum/Maximum number of cells matated in one feature, optional (default: None)
    - min.max_pct_cells (float): Minimum/Maximum percentage of cells mutated in one feature, optional (default: None)
        
*return:  None 
```
<a name="filter_cells"/>

> 10. Filter out cells based on different metrics.
```
filter_cells(adata,min_n_features=None, max_n_features=None, min_pct_features=None, max_pct_features=None):
    
 * Parameters:
    - adata           (AnnData): anndata object that contains the mutation / gene expression profiles.
    - min/max_n_cells     (int): Minimum/Maximum number of mutations in the cell, optional (default: None)
    - min.max_pct_cells (float): Minimum/Maximum percentage of mutations in the cell, optional (default: None)
        
*return:  None 
```




