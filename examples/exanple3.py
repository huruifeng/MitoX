import mitox as mx

import matplotlib
matplotlib.use('agg')

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
np.random.seed(12)

bam_folder = "example/E-MTAB-5061/bam"
mut_adata = mx.gen_mut_profile(bam_folder,sp="human",seq_type="sc",combined_bam=False)
mx.add_metadata(mut_adata,file_name="example/E-MTAB-5061/metadata.txt")
mx.write_table(mut_adata,file="T2Dmut.txt",sep="\t")

mx.filter_features(mut_adata, min_n_cells=10)
mx.filter_cells(mut_adata, min_n_features=10)
mx.plot_mean_std(mut_adata)

mx.select_top_varible_variants(mut_adata,top_n = 50)
mx.plot_clustermap(mut_adata,features="var_muts",ann_color=["individual","cell type","disease"], ann_label=None)

mx.select_top_principal_components(mut_adata,features="var_muts",n_pc=20, fig_name="var_muts_pc.png")
mx.select_top_principal_components(mut_adata)

mx.dimension_reduction(mut_adata)
mx.plot_clusters(mut_adata)
mx.plot_clusters(mut_adata,fig_name="cluster_tsne.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_disease.png", cluster_color="disease")

mx.dimension_reduction(mut_adata,dim_reducer="PCA")
mx.plot_clusters(mut_adata,fig_name="cluster_pca.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_pca_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_pca_disease.png", cluster_color="disease")

mx.dimension_reduction(mut_adata,dim_reducer="UMAP")
mx.plot_clusters(mut_adata,fig_name="cluster_umap.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_umap_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_umap_disease.png", cluster_color="disease")

mx.dimension_reduction(mut_adata,features="var_muts")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_var_muts.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_var_muts_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_var_muts_disease.png", cluster_color="disease")

mx.dimension_reduction(mut_adata,features="top_pcs")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_top_pcs.png", cluster_color="individual")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_top_pcs_cell.png", cluster_color="cell type")
mx.plot_clusters(mut_adata,fig_name="cluster_tsne_top_pcs_disease.png", cluster_color="disease")

# bam_folder = "example/example/10x_1"
expr_adata = mx.gen_expression_profile(bam_folder,sp="human",gtf_file="example/E-MTAB-5061/Homo_sapiens.GRCh38.101.gtf.gz")
expr_df = expr_adata.to_df().T

###
mut_adata = mx.read_table("mut_t2d.txt",sep="\t",type="mut")
expr_adata = mx.read_table("expr_t2d.txt",sep="\t",type="expr")
mx.add_metadata(expr_adata,file_name="example/E-MTAB-5061/metadata.txt")

combined_adata = mx.combine_mut_expr(mut_adata,expr_adata)

mx.dimension_reduction(combined_adata)
mx.plot_clusters(combined_adata,fig_name="cluster_tsne_gene_MT-ATP8.png", cluster_color="MT-ATP8")

mx.plot_expr(expr_adata,gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)
mx.plot_expr(expr_adata,group="label",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)
mx.plot_expr(expr_adata,group="disease",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)
mx.plot_expr(expr_adata,group="disease",using={"cell type":["delta cell"]},gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)

mx.plot_expr(expr_adata,gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False)
mx.plot_expr(expr_adata,group="label",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False)
mx.plot_expr(expr_adata,group="disease",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False)
mx.plot_expr(expr_adata,group="disease",using={"cell type":["delta cell"]},gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True,all_in_one=False)

mx.plot_expr(combined_adata,gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)
mx.plot_expr(combined_adata,group="label",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)
mx.plot_expr(combined_adata,group="disease",gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)
mx.plot_expr(combined_adata,group="disease",using={"cell type":["delta cell"]},gene=["MT-ATP6","MT-ATP8","MT-CO1","MT-CO2"],plot_type="violin",strip=True)

