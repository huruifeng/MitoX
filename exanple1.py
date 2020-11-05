import mitox as mx

import numpy as np
np.random.seed(12)

bam_folder = "example/GSE74310/bam"
samples = mx.read_sam(bam_folder,sp="human",seq_type="sc",combined_bam=False)

adata = mx.gen_mut_profile(samples)
mx.add_metadata(adata,"example/GSE74310/metadata.txt")

mx.filter_features(adata, min_n_cells=10)
mx.filter_cells(adata, min_n_features=10)
mx.plot_mean_std(adata)

mx.select_top_varible_variants(adata,top_n = 50)

mx.cal_distance(adata)

mx.plot_clustermap(adata,features="var_muts",ann_color=["sample","cell_type"], ann_label=None)
mx.plot_clustermap(adata,features="distance",ann_color=["sample","cell_type"], ann_label=None, fig_name="cluster_map_dis.png")

mx.select_top_principal_components(adata)
mx.select_top_principal_components(adata,features="var_muts",n_pc=20, fig_name="var_muts_pc.png")

mx.dimension_reduction(adata)
mx.plot_clusters(adata)

mx.dimension_reduction(adata,dim_reducer="PCA")
mx.plot_clusters(adata,fig_name="cluster_pca.png", cluster_color="label")

mx.dimension_reduction(adata,dim_reducer="UMAP")
mx.plot_clusters(adata,fig_name="umap_cluster.png",cluster_n=4)
mx.plot_clusters(adata,fig_name="umap_cluster_label.png", cluster_color="label")
mx.plot_clusters(adata,fig_name="umap_cluster_sample.png", cluster_color="sample")
mx.plot_clusters(adata,fig_name="umap_cluster_cell.png", cluster_color="cell_type")
mx.plot_clusters(adata,fig_name="umap_cluster_n_muts.png", cluster_color="n_muts")


mx.dimension_reduction(adata,n_components=3, dim_reducer="UMAP")
mx.plot_clusters(adata,fig_name="umap_cluster_3d.png",cluster_n=4)
mx.plot_clusters(adata,fig_name="umap_cluster_label_3d.png", cluster_color="label")
mx.plot_clusters(adata,fig_name="umap_cluster_sample_3d.png", cluster_color="sample")
mx.plot_clusters(adata,fig_name="umap_cluster_cell_3d.png", cluster_color="cell_type")
mx.plot_clusters(adata,fig_name="umap_cluster_n_muts_3d.png",cluster_color="n_muts",color_pal="flare")





