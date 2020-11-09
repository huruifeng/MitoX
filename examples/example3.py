import mitox as mx
import pandas as pd

bam_folder = "example3/bam"
samples = mx.read_sam(bam_folder,seq_type="sc",combined_bam=True,min_read=100)
x_data = mx.gen_mut_profile(samples)

## using all variants, It takes time and memery.
# mx.plot_clustermap(x_data,features="all",fig_name="all_mut_clustermap.png")

mut_mean_sd = mx.mut_mean_std(x_data.to_df().T)
top_mean_std = mut_mean_sd.nlargest(20,"std",keep="all")
top_mean_std.to_csv("top20_std_variants.txt",sep="\t")

mx.plot_mean_std(x_data,fig_name="all_mut_mean_sd.png")
mx.select_top_varible_variants(x_data,top_n=20, fig_name="top20_mean_sd.png")

mx.plot_clustermap(x_data,fig_name="top20Mut_clustermap.png")

mx.add_metadata(x_data,file_name="example3/metadata.txt")

## after adding metadata, we can annotate the columns of the clustermap
mx.plot_clustermap(x_data,ann_color=["label"], ann_label=["Cell type"],fig_name="top20Mut_clustermap_cell_anno.png")

mx.dimension_reduction(x_data)
## if  cluster_color is not set, the cluster will be detcted using algorithm DBSCAN.
## User can also set this algrorithm to KMeans, ane set the number of clusters
mx.plot_clusters(x_data,fig_name="cluster_default.png")
mx.plot_clusters(x_data,fig_name="cluster_tsne.png", cluster_color="label")

##PCA
mx.dimension_reduction(x_data,dim_reducer="PCA")
mx.plot_clusters(x_data,fig_name="cluster_pca.png", cluster_color="label")