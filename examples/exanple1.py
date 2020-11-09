## import MitoX and necessary python packages
import mitox as mx

import matplotlib
matplotlib.use('agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
np.random.seed(12)

## Read bulk ATAC-seq bam files for C9, D6, G10 cell clones
bam_folder = "example1/C9D6G10/bams/ATAC"
ATAC_samples = mx.read_sam(bam_folder,sp="human",seq_type="bulk",combined_bam=False)
## Generate the heteroplasmy profiles
atac_adata = mx.gen_mut_profile(ATAC_samples)
## Write heteroplasmy profiles to local files
mx.write_table(atac_adata,file="ATACmut.txt",sep="\t")
## Transfer AnnData format to pandas DataFrame format
atac_df = atac_adata.to_df().T
## Calculate the coverage in G10 clone bulk ATAC-seq data
G10_atac_coverage = mx.get_coverage(atac_adata,name="G10_ATAC")

## Read bulk RNA-seq bam files for C9, D6, G10 cell clones
bam_folder = "example1/C9D6G10/bams/RNA"
RNA_samples = mx.read_sam(bam_folder,sp="human",seq_type="bulk",combined_bam=False)
rna_adata = mx.gen_mut_profile(RNA_samples)
mx.write_table(rna_adata,file="RNAmut.txt",sep="\t")
rna_df = rna_adata.to_df().T
## Calculate the coverage in G10 clone bulk RNA-seq data
G10_rna_coverage = mx.get_coverage(rna_adata,name="G10_RNA")

## Combine the G10 ATAC-seq and RNA-seq coverage information into a pd.DataFrame for plot
G10_cov_df = pd.concat([G10_atac_coverage,G10_rna_coverage], axis=1)
## Plot the coverage information in a circle plot(default)
mx.plot_coverage(G10_cov_df)
## Plot the coverage information in a plane plot(default)
mx.plot_coverage(G10_cov_df,fig_name="coverage_plane.png",plot_type="plane")

##Filter variants, only keep those with heteroplasmies > 0.01
C9_df = pd.DataFrame([atac_df["C9_ATAC"],rna_df["C9_RNA"]]).T
C9_df = C9_df[(C9_df>0.01).all(axis=1)]

D6_df = pd.DataFrame([atac_df["D6_ATAC"],rna_df["D6_RNA"]]).T
D6_df = D6_df[(D6_df>0.01).all(axis=1)]

G10_df = pd.DataFrame([atac_df["G10_ATAC"],rna_df["G10_RNA"]]).T
G10_df = G10_df[(G10_df>0.01).all(axis=1)]

##Plot the heteroplasmies between ATAC-seq data and RNA-seq data
fig = plt.figure(figsize=(4,4))
ax = fig.add_subplot()

ax.scatter(C9_df["C9_ATAC"],C9_df["C9_RNA"],c="red", label="C9", s=10)
ax.scatter(D6_df["D6_ATAC"],D6_df["D6_RNA"],c="green",label="D6",s=10)
ax.scatter(G10_df["G10_ATAC"],G10_df["G10_RNA"],c="blue",label="G10",s=10)
ax.plot([0,1],[0,1], dashes=[6, 3], c=".3")

ax.set_xlim(-0.05, 1.05)
ax.set_ylim(-0.05, 1.05)
ax.set_xlabel("Allel frequency - ATAC")
ax.set_ylabel("Allel frequency - RNA")
ax.legend(loc="lower right",title="cell clone",frameon=False)
plt.tight_layout()
fig.savefig("ATAC_RNA.png")
plt.close(fig)

##############
## Processing single cell sequencing data
# Read BAM files, each bam file is for one cell
bam_folder = "example1/C9D6G10/bams/scATAC"
scATAC_samples = mx.read_sam(bam_folder,sp="human",seq_type="sc",combined_bam=False)
scatac_adata = mx.gen_mut_profile(scATAC_samples)
mx.add_metadata(scatac_adata,file_name="example1/C9D6G10/scATAC_meta.txt")
mx.write_table(scatac_adata,file="scATACmut.txt",sep="\t")
scatac_df = scatac_adata.to_df().T

bam_folder = "example1/C9D6G10/bams/scRNA"
scRNA_samples = mx.read_sam(bam_folder,sp="human",seq_type="sc",combined_bam=False)
scrna_adata = mx.gen_mut_profile(scRNA_samples)
mx.add_metadata(scrna_adata,file_name="example1/C9D6G10/scRNA_meta.txt")
mx.write_table(scrna_adata,file="scRNAmut.txt",sep="\t")
scrna_df = scrna_adata.to_df().T

## selected clone specific variants
sc_atac_D6_ls = mx.select_specific_variants(scrna_adata,sample_name="D6")
sc_atac_C9_ls = mx.select_specific_variants(scrna_adata,sample_name="C9")
sc_atac_G10_ls = mx.select_specific_variants(scrna_adata,sample_name="G10")
sc_atac_var_ls = sc_atac_D6_ls+sc_atac_C9_ls+sc_atac_G10_ls

sc_rna_D6_ls = mx.select_specific_variants(scrna_adata,sample_name="D6")
sc_rna_C9_ls = mx.select_specific_variants(scrna_adata,sample_name="C9")
sc_rna_G10_ls = mx.select_specific_variants(scrna_adata,sample_name="G10")
sc_rna_var_ls = sc_rna_D6_ls+sc_rna_C9_ls+sc_rna_G10_ls

union_ls = list(set(sc_atac_var_ls+sc_atac_var_ls))

# Save to local files
var_atac_df = scatac_df.loc[union_ls,:]
var_atac_df.to_csv("clone_specific_scATACseq.txt",sep="\t")
var_rna_df = scrna_df.loc[union_ls,:]
var_atac_df.to_csv("clone_specific_scRNAseq.txt",sep="\t")

## transfer to percentage
var_atac_df = var_atac_df * 100
var_rna_df = var_rna_df * 100

## Plot and compare the clone specific variants between scATAC-seq and scRNA-seq in the three cell clones
fig = plt.figure(num=1, figsize=(20, 8),dpi=80)
i=0
for var_i in union_ls:
    i+=1
    ax = fig.add_subplot(2, 5, i)
    if var_i in sc_atac_C9_ls:
        bulk_atac = atac_df.loc[var_i, "C9_ATAC"]
        bulk_atac_str = "ATAC:" + str(round(bulk_atac*100, 2))
    if var_i in sc_atac_D6_ls:
        bulk_atac = atac_df.loc[var_i,"D6_ATAC"]
        bulk_atac_str = "ATAC:"+str(round(bulk_atac*100, 2))
    if var_i in sc_atac_G10_ls:
        bulk_atac = atac_df.loc[var_i,"G10_ATAC"]
        bulk_atac_str = "ATAC:"+str(round(bulk_atac*100, 2))

    if var_i in sc_rna_C9_ls:
        bulk_rna = rna_df.loc[var_i,"C9_RNA"]
        bulk_rna_str = "RNA:"+str(round(bulk_rna*100, 2))
    if var_i in sc_rna_D6_ls:
        bulk_rna = rna_df.loc[var_i,"D6_RNA"]
        bulk_rna_str = "RNA:"+str(round(bulk_rna*100, 2))
    if var_i in sc_rna_G10_ls:
        bulk_rna = rna_df.loc[var_i,"G10_RNA"]
        bulk_rna_str = "RNA:"+str(round(bulk_rna*100, 2))

    scrna_values_C9 = var_rna_df.loc[var_i, var_rna_df.columns.str.contains("C9")]
    scrna_values_D6=var_rna_df.loc[var_i,var_rna_df.columns.str.contains("D6")]
    scrna_values_G10 = var_rna_df.loc[var_i, var_rna_df.columns.str.contains("G10")]

    scatac_values_C9 = var_atac_df.loc[var_i, var_atac_df.columns.str.contains("C9")]
    scatac_values_D6 = var_atac_df.loc[var_i, var_atac_df.columns.str.contains("D6")]
    scatac_values_G10 = var_atac_df.loc[var_i, var_atac_df.columns.str.contains("G10")]

    data = {"Assay":[],"Cell":[],"Value":[]}
    data["Value"] = list(scrna_values_C9)+list(scrna_values_D6)+list(scrna_values_G10) + \
                    list(scatac_values_C9)+list(scatac_values_D6)+list(scatac_values_G10)
    data["Cell"] = ["C9"] * len(scrna_values_C9) + ["D6"] * len(scrna_values_D6) + ["G10"] * len(scrna_values_G10) + \
                   ["C9"] * len(scatac_values_C9) + ["D6"] * len(scatac_values_D6) + ["G10"] * len(scatac_values_G10)
    data["Assay"] = ["scRNA"] * (len(scrna_values_C9) + len(scrna_values_D6) + len(scrna_values_G10)) + \
                    ["scATAC"] *(len(scatac_values_C9) + len(scatac_values_D6) + len(scatac_values_G10))

    data_df = pd.DataFrame(data=data)
    sns.set_style("ticks")
    sns.boxplot(x="Assay",y="Value",hue="Cell",data=data_df,ax=ax)
    sns.stripplot(x="Assay",y="Value",hue="Cell", data=data_df,ax=ax,dodge=True)
    ax.get_legend().remove()
    ax.set_title(var_i + " ("+bulk_rna_str+" | "+bulk_atac_str+")")
    ax.set_xlabel("")
    ax.set_ylabel("")

fig.text(0, 0.5, '% heteroplasmy', ha='center', va='center', rotation='vertical')
handles, labels = ax.get_legend_handles_labels()
n = len(data_df["Cell"].unique())
l = plt.legend(handles[0:n], labels[0:n], bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.,frameon=False)

fig.tight_layout()  # 调整整体空白
plt.subplots_adjust(wspace=0.2, hspace=0.2)  # 调整子图间距

fig.savefig("boxplot_clone_specific_variants.png")
plt.close(fig)


