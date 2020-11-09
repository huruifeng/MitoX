import mitox as mx

import matplotlib
matplotlib.use('agg')

import numpy as np
np.random.seed(12)

# Read BAM files, each bam file is for one cell
bam_folder = "example2/E-MTAB-5061/bam"
mut_adata = mx.gen_mut_profile(bam_folder,sp="human",seq_type="sc",combined_bam=False,n_jobs=4)
mx.add_metadata(mut_adata,file_name="example2/E-MTAB-5061/metadata.txt")
mx.write_table(mut_adata,file="T2Dmut.txt",sep="\t")

## Filter cells and variants
mx.filter_features(mut_adata, min_n_cells=10)
mx.filter_cells(mut_adata, min_n_features=10)
# plot the distribution of Standard Deviation and Mean values of all heteroplasmies
mx.plot_mean_std(mut_adata)

## select top_varible_variants, and plot a clustermap using these variants
mx.select_top_varible_variants(mut_adata,top_n = 50)
mx.plot_clustermap(mut_adata,features="var_muts",ann_color=["individual","cell type","disease"], ann_label=None)

## select and plot top_principal_components
mx.select_top_principal_components(mut_adata,features="var_muts",n_pc=20, fig_name="var_muts_pc.png")
mx.select_top_principal_components(mut_adata)

## Dimension reduction using all variant information
##tSNE(default)
mx.dimension_reduction(mut_adata)
## if  cluster_color is not set, the cluster will be detcted using algorithm DBSCAN.
## User can also set this algrorithm to KMeans, ane set the number of clusters
mx.plot_clusters(mut_adata,fig_name="cluster_default.png")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_disease.png", cluster_color="disease")

##PCA
mx.dimension_reduction(mut_adata,dim_reducer="PCA")
mx.plot_clusters(mut_adata,fig_name="cluster_pca.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_pca_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_pca_disease.png", cluster_color="disease")

## UMAP
mx.dimension_reduction(mut_adata,dim_reducer="UMAP")
mx.plot_clusters(mut_adata,fig_name="cluster_umap.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_umap_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_umap_disease.png", cluster_color="disease")

## Dimension reduction using top variable variants
mx.dimension_reduction(mut_adata,features="var_muts")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_var_muts.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_var_muts_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_var_muts_disease.png", cluster_color="disease")

## Dimension reduction using top principal components
mx.dimension_reduction(mut_adata,features="top_pcs")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_top_pcs.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_top_pcs_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_top_pcs_disease.png", cluster_color="disease")

## Generate gene expression profiles
expr_adata = mx.gen_expression_profile(bam_folder,sp="human",gtf_file="example2/E-MTAB-5061/Homo_sapiens.GRCh38.101.gtf.gz",n_jobs=4)
expr_df = expr_adata.to_df().T
expr_df.to_csv("T2Dexpr.txt",sep="\t")

### Read data from table
# mut_adata = mx.read_table("T2Dmut.txt",sep="\t",type="mut")
# expr_adata = mx.read_table("T2Dexpr.txt",sep="\t",type="expr")

mx.plot_expr(expr_adata,gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True, fig_name="expr_plot_in_one.png")
mx.plot_expr(expr_adata,gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False,fig_name="expr_plot_separate.png")

mx.add_metadata(expr_adata, file_name="example2/E-MTAB-5061/metadata.txt")
mx.plot_expr(expr_adata,group="label",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,fig_name="expr_plot_label.png")
mx.plot_expr(expr_adata,group="disease",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,fig_name="expr_plot_disease.png")
mx.plot_expr(expr_adata,group="disease",using={"cell type":["delta cell"]},gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,fig_name="expr_plot_disease_cell.png")

mx.plot_expr(expr_adata,group="label",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False,fig_name="expr_plot_label_separate.png")
mx.plot_expr(expr_adata,group="disease",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False,fig_name="expr_plot_disease_separate.png")
mx.plot_expr(expr_adata,group="disease",using={"cell type":["delta cell"]},gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False,fig_name="expr_plot_disease_separate_cell.png")

## Using the combined anndata to plot gene expression patterns
combined_adata = mx.combine_mut_expr(mut_adata,expr_adata)

mx.dimension_reduction(combined_adata)
mx.plot_clusters(combined_adata,fig_name="cluster_tsne_gene_MT-ATP8_combine.png", cluster_color="MT-ATP8")

mx.plot_expr(combined_adata,gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,fig_name="expr_combine.png")
mx.plot_expr(combined_adata,group="label",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,fig_name="expr_combine_label.png")
mx.plot_expr(combined_adata,group="disease",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,fig_name="expr_combine_disease.png")
mx.plot_expr(combined_adata,group="disease",using={"cell type":["delta cell"]},gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,fig_name="expr_combine_disease_cell.png")

## tSNE plot on expresion data
import pandas as pd
from sklearn.manifold import TSNE
import seaborn as sns
import matplotlib.pyplot as plt

expr_rpkm = pd.read_csv("example2/E-MTAB-5061/rpkm.txt",sep="\t",index_col=0, header=0)
expr_counts = pd.read_csv("example2/E-MTAB-5061/counts.txt",sep="\t",index_col=0,header=0)

expr_data = expr_rpkm

meta_data = pd.read_csv("example2/E-MTAB-5061/metadata.txt",sep="\t",index_col=0, header=0)

meta_data = meta_data.loc[meta_data['disease']=="normal",:]

samples = meta_data.index
expr_data = expr_data.loc[:,samples]

tsne = TSNE(n_components=2, perplexity=30.0, n_iter=1000,random_state=12)
X_embedded = tsne.fit_transform(expr_data.T)

sns.set(style="ticks")
f, ax = plt.subplots(figsize=(12, 9))
sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=meta_data["cell type"],legend="full", palette="husl",ax =ax)
ax.set(xlabel="tSNE 1")
ax.set(ylabel="tSNE 2")
ax.legend(loc=2, bbox_to_anchor=(1, 1), frameon=False)

plt.tight_layout()
f.savefig("all_gene_expr_tsne_normal_cell.png")
plt.close()


