
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
    * [plot_mean_std()](#plot_mean_std)
    * [cal_distance()](#cal_distance)
    * [write_distance()](#write_distance)
    * [select_specific_variants()](#select_specific_variants)
    * [select_top_varible_variants()](#select_top_varible_variants)
    * [select_top_principal_components()](#select_top_principal_components)
    * [plot_clustermap()](#plot_clustermap)
    * [dimension_reduction()](#dimension_reduction)
    * [combine_mut_expr()](#combine_mut_expr)
    * [plot_expr()](#plot_expr)


    

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
```python
 mut_adata = gen_mut_profile(bams, fmt="bam", sp="human", chr_name="MT", seq_type="sc", combined=False, tag_name="CB",
                           barcodes="NULL", min_read=200, max_depth=1e5, min_baseq=25, min_mapq=0, n_jobs=2)
    
 * Parameters:
    - bams      (str): Path pointing to BAM alignment file(s).
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
```python
 expr_adata = gen_expression_profile(bams, gtf_file, fmt="bam", sp="human", chr_name="MT", eq_type="sc",
                                    combined=False, tag_name="CB", barcodes="NULL", min_read=200, stranded=True,
                                    feature_type="exon", min_mapq=10, n_jobs=2):
    
 * Parameters:
    - bams      (str): Path pointing to BAM alignment file(s).
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
```python
 samples = read_sam(bam_folder, fmt="bam",sp="human", chr_name="MT", seq_type="sc",combined=False, tag_name="CB", barcodes="NULL",
                       min_read=200, min_baseq=25, min_mapq=0, n_jobs=2):
    
 * Parameters:
    - bam_folder(str): Path pointing to BAM alignment file(s).
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
*return (Sample): A sample or sample list that contiains coverage information for each sample.
```

<a name="read_table"/>

> 4. Read data from text file.
```python
 x_adata = read_table(file,sep="\t",type="mut"):
    
 * Parameters:
    - file  (str): Path pointing to file.
    - sep   (str): Delimiter to use. 
    - type  (str): Data type - mutation profile table ('mut') or gene expression profile table ('expr').
*return (AnnData): Anndata object that contains the mutation / gene expression profiles
```
<a name="write_table"/>

> 5. Write AnnData into local file.
```python
write_table(adata,file,sep="\t"):
    
 * Parameters:
    - adata (AnnData): Anndata object that contains the mutation / gene expression profiles.
    - file      (str): File path and name 
    - sep       (str): String of length 1. Field delimiter for the output file..
*return : None
```

<a name="get_coverage"/>

> 6. Get the coverage information for sample(s).
```python
cov = get_coverage(adata,name="all")
    
 * Parameters:
    - adata (AnnData): Anndata object that contains the mutation / gene expression profiles.
    - name      (str): Sample name, when "all", returns a dataframe for all samples 
*return (Series or DataFrame): coverage data for sample(s) 
```
<a name="plot_coverage"/>

> 7. Plot the coverage information for sample(s).
```python
plot_coverage(x_df,color=None,fig_name="coverage_plot.png",fig_size=None,plot_type="circle",log2=True)
    
 * Parameters:
    - x_df (Series,DataFrame): A Series of DataFrame that contains the coverage data. 
                               It can be the return value from get_coverage().
    - color       (list, str): A color name (e.g. 'red') for the Series，
                               OR, a list of colors for the DataFrome.
    - fig_name          (str): Figure loaction and name for saving.
    - fig_size (float, float): Width, height in inches.
    - plot_type         (str): "circle" or "plane".
    - log2             (bool): Do log2 transform on coverage values or not.
*return:  None 
```
<a name="add_metadata"/>

> 8. Add annotation data (metadata) to samples.
```python
add_metadata(adata,file_name,delimiter='\t')
    
 * Parameters:
    - adata  (AnnData): Anndata object that contains the mutation / gene expression profiles.
    - file_name  (str): File path and name to the metadata file. 
    - sep        (str): String of length 1. Field delimiter of the metadata file.
*return:  None 
```
<a name="filter_features"/>

> 9. Filter out features based on different metrics.
```python
filter_features(adata, min_n_cells = None, max_n_cells=None, min_pct_cells = None, max_pct_cells=None)
    
 * Parameters:
    - adata           (AnnData): Anndata object that contains the mutation / gene expression profiles.
    - min/max_n_cells     (int): Minimum/Maximum number of cells matated in one feature, optional (default: None)
    - min.max_pct_cells (float): Minimum/Maximum percentage of cells mutated in one feature, optional (default: None)     
*return:  None 
```
<a name="filter_cells"/>

> 10. Filter out cells based on different metrics.
```python
filter_cells(adata,min_n_features=None, max_n_features=None, min_pct_features=None, max_pct_features=None):
    
 * Parameters:
    - adata           (AnnData): anndata object that contains the mutation / gene expression profiles.
    - min/max_n_cells     (int): Minimum/Maximum number of mutations in the cell, optional (default: None)
    - min.max_pct_cells (float): Minimum/Maximum percentage of mutations in the cell, optional (default: None)      
*return:  None 
```
<a name="plot_mean_std"/>

> 11. Plot the mean, standard deviation distribution.
```python
 plot_mean_std(adata, fig_name = "mean_std.png")
    
 * Parameters:
    - adata     (AnnData): Anndata object that contains the mutation / gene expression profiles.
    - fig_name      (str): Figure loaction and name for saving.       
*return:  None 
```
<a name="cal_distance"/>

> 12. Calculate the relatedbess of sample using heteroplasmy profile.

The distance is defined as: <img src="https://latex.codecogs.com/gif.latex?d_{ij}=\frac{\sum_{x}^{}(\sqrt{\left&space;|&space;AF_{x,j=i}-&space;AF_{x,j}\right&space;|}*(1_{c_{x,i}>0}*1_{c_{x,j}>0}))}{\sum_{x}^{}(1_{c_{x,i}>0}*1_{c_{x,j}>0})}" />
<!--img src="http://chart.googleapis.com/chart?cht=tx&chl=d_{ij}=\frac{\sum_{x}^{}(\sqrt{|AF_{x,i}AF_{x,j}|}*(1_{c_{x,i}>0}*1_{c_{x,j}>0}))}{\sum_{x}^{}(1_{c_{x,i}>0}*1_{c_{x,j}>0})}" style="border:none;"-->

```python
 cal_distance(adata,features="all")
    
 * Parameters:
    - adata (AnnData): Anndata object that contains the heteroplasmy profile.
    - features  (str): Feature that will be used fro the calcualtion,'var_muts' or 'all'.   
*return:  None 
```

<a name="write_distance"/>

> 13. Write relatedbess matrix to local file.
```python
write_distance(adata, file,sep="\t")
    
 * Parameters:
    - adata (AnnData): Anndata object that contains the heteroplasmy profile.
    - file      (str): File path and name to save.
    - sep       (str): String of length 1. Field delimiter for the output file.

*return : None
```

<a name="select_specific_variants"/>

> 14. Select clone or sample specific variants.
```python
select_specific_variants(adata,sample_name, min_af = 0.01, percent=0.8, check_other=True,other_min_af=0.01,other_n=5, other_percent=0)
    
 * Parameters:
    - adata    (AnnData): Anndata object that contains the heteroplasmy profile.
    - sample_name  (str): Sample name that will be used to extract the specific variants.
    - min_af     (float): The minimun allele frequency(AF) that will be used to filter the variants.
    - percent    (float): The percent of sample that have this variants with allele frequency(AF) > min_af.
    - check_other (bool): Apply the filter conditions to the other samples. if False, the result may contain homoplamic variants.
    - other_min_af, other_n ,other_percent: Used to define the conditons to filter variants other than the selected samples.
                                            other_n or other_percent will be used whichever defines the smaller sample number. 
*return (list): The list variant IDs.
```
<a name="select_top_varible_variants"/>

> 15. Show top varible variants in a scatter plot, also store the information in tha adata object.
```python
select_top_varible_variants(adata, top_n = 50, percentile=20, fig_name = 'std_vs_means.png',fig_size = (8, 8),pad = 1.08)    
 
* Parameters:
    - adata    (AnnData): Anndata object that contains the heteroplasmy profile.
    - top_n        (int): The number of variants will be selected.
    - percentile   (int): The top percentile of variants will be selected.
    - fig_name     (str): Figure name for saving.
    - fig_size (float, float): Width, height in inches.
    - pad        (float): Padding between the frame and figure edge.  
*return: None.
```
<a name="select_top_principal_components"/>

> 16. Select top principal components and store the information in tha adata object for further use.
```python
select_top_principal_components(adata, features="var_muts", n_pc=50, max_pc=100, first_pc=True,fig_name='top_pcs.png',fig_size=(8, 8),pad=1.08)  
 
* Parameters:
    - adata    (AnnData): Anndata object that contains the heteroplasmy profile.
    - features     (str): To define the dataset that will be used for principal components analysis.
    - n_pc         (int): The top n PCs will be selected and stored in adata.
    - first_pc    (bool): If the first PC will be used ot not.
    - fig_name     (str): Figure name for saving.
    - fig_size (float, float): Width, height in inches.
    - pad        (float): Padding between the frame and figure edge.  
*return: None.
```
<a name="plot_clustermap"/>

> 17. Select top principal components and store the information in tha adata object for further use.
```python
plot_clustermap(adata, features="var_muts", ann_color=None, ann_label=None, fig_name="clusterheatmap.png", vmin=0, vmax=1,
                fig_size=(12,12), xlabel=False,ylabel=True,xlabel_size=8,ylabel_size=8,xcluster=True, ycluster=True)
 
* Parameters:
    - adata    (AnnData): Anndata object that contains the heteroplasmy profile.
    - features     (str): To define the dataset that will be used for ploting the clustermap,"all" or "mut_vars".
    - ann_color   (list): annotation labels that will be used the columns color annotation.
    - ann_label   (list): Rename the values in ann_color.
    - fig_name     (str): Figure name for saving.
    - fig_size (float, float): Width, height in inches.
    - vmin/vmax             (float): 
    - xlabele/ylabel         (bool): If show xlabele/ylabel or not.
    - xlabel_size/ylabel_size (int): Set xlabel font size/ylabel font size
    - xcluster/ycluster      (bool): If do cluster for columns and rows or not.
*return: None.
```
<a name="dimension_reduction"/>

> 18. Do the dimension reduction on samples.
```python
dimension_reduction(adata,features="all", n_neighbors=20, n_components=2, dim_reducer="tsne")
 
* Parameters:
    - adata    (AnnData): Anndata object that contains the heteroplasmy profile.
    - features     (str): To define the dataset that will be used for dimension reduction, "all", "mut_vars","top_pcs", or "distance".
    - n_neighbors  (int): Used for UMAP() method. 
    - n_components (int): Number of components to keep. Typically, it is 2, or 3. 
    - dim_reducer  (str): Method used to do the dimension reduction, "tSNA","PCA",or "UMAP".
*return: None.
```
<a name="plot_clusters"/>

> 19. Plot clusters on the results of dimension_reduction().
```python
plot_clusters(adata, fig_name="cluster.png", cluster_color="auto",method="DBSCAN",cluster_n=2,color_pal="Blues",cbar_frac=0.025, log2=False)
 
* Parameters:
    - adata    (AnnData): Anndata object that contains the heteroplasmy profile.
    - fig_name     (str): Figure name for saving.
    - cluster_color(str): Set cluster colors, if "auto" the cluster will be detected using "method" (DBSCAN, kmeans), 
                          Otherwise, the value should be one of the column names in metadata.
    - method       (str): Method used to detect if cluster:
    - cluster_n    (int): When method="kmeans", set the number of cluster that will be detected.
    - color_pal    (int): When cluster_color="[gene_name]", set the color palette for continuous data.
    - cbar_frac  (float): Color bar size.
    - log2        (bool): When cluster_color="[gene_name]", do log2 transformation.
*return: None.
```
<a name="combine_mut_expr"/>

> 19. Combine the heteroplasmy profile with gene exprssion profile.
```python
combine_mut_expr(mut_adata, expr_adata) 

* Parameters:
    - mut_adata  (AnnData): Anndata object that contains the heteroplasmy profile.
    - expr_adata (AnnData): Anndata object that contains the gene expression profile..
*return (AnnData): Combined anndata object.
```
<a name="plot_expr"/>

> 20. Combine the heteroplasmy profile with gene exprssion profile.
```python
plot_expr(adata, gene, group=None, using={}, fig_name="expr_plot.png", fig_size=(12, 8), plot_type="boxplot",
          strip=False,log2=True, all_in_one=True,col_n=0,xrotate=30, **kwargs)

* Parameters:
    - adata   (AnnData): Anndata object that contains the heteroplasmy profile.
    - gene       (list): A list of genes to plot.
    - group       (str): Define the groups to compare, it should be one of the column names in metadata.
    - using      (dict): A dict to define which the sub-dataset that will be ploted.
    - fig_name    (str): Figure name for saving.
    - fig_size (float, float): Width, height in inches.
    - strip      (bool): Plot the stripplot.
    - log2       (bool): Do log2 transformation of values.
    - all_in_one (bool): Plot all genes in one integrated plot or draw each gene in a subplot.
    - col_n       (int): If all_in_one is False, set the number of subplot in each row.
    - xrotate     (int): Set the rotsted degree of x labels.
    - **kwargs         : Other options for sns.boxplot(), sns.violinplot(), sns.stripplot(), such as "color='red'".
*return: None.
```














