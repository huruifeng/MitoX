from .bamreader import read_sam
from .global_data import human_mt_seq
from .global_data import human_mt_len
from .global_data import mus_mt_seq
from .global_data import mus_mt_len

from .sample import Sample
from .funcs import *

import os
import pandas as pd
import numpy as np
import anndata as ad

from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import umap

from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans


import matplotlib
matplotlib.use('agg')

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

np.random.seed(12)

def read_table(file,sep="\t",type="mut"):
    df = pd.read_csv(file,sep=sep,index_col=0,header=0)
    adata = ad.AnnData(df.T)
    adata.uns[type] = True
    return adata

def write_table(adata,file,sep="\t"):
    x_df = adata.to_df().T
    x_df.to_csv(file,sep=sep)
    return None

def get_coverage(adata,name="all"):
    cov_df = adata.varm["coverage"]
    index_ls = cov_df.index.str[:-3]
    cov_df["pos_i"] = index_ls
    cov_df.drop_duplicates(keep="first", inplace=True)
    cov_df.drop(['pos_i'], axis = 1, inplace = True)
    cov_df.index = cov_df.index.str[:-3]
    if name == "all":
        return cov_df
    else:
        if name not in cov_df.columns:
            print("ERROR: There is no sample name: "+name)
        else:
            return cov_df[name]


def plot_coverage(x_df,color=None,fig_name="coverage_plot.png",fig_size=None,plot_type="circle",log2=True):
    if isinstance(x_df,pd.Series):
        x_df = pd.DataFrame(x_df)

    cov_df = x_df
    name = cov_df.columns
    pos_len = x_df.shape[0]
    if color==None or (isinstance(color,list) and len(color) != len(name)):
        color = sns.color_palette("hls", n_colors=len(name)).as_hex()

    pos_ls = np.array(cov_df.index, dtype=np.float)
    pos_len_quater = int(pos_len / 4)
    pos_label = ['0', str(pos_len_quater),str(pos_len_quater * 2),str(pos_len_quater * 3)]
    ylabel_str = "Coverage"
    if log2:
        cov_df = np.log2(cov_df + 1)
        ylabel_str="log2(coverage + 1)"


    if fig_size == None:
        fig_size = (8,8) if plot_type == "circle" else (12,3)

    fig = plt.figure(figsize=fig_size)
    if plot_type == "circle":
        pos_ls = pos_ls / pos_len
        theta = 2 * np.pi * pos_ls

        ax = fig.add_subplot(111, polar=True)
        for i,name_i in enumerate(name):
            c_i = color[i]
            r = cov_df[name_i]
            ax.plot(theta, r, c=c_i,label=name_i)

        # ax.set_xticklabels(pos_label)
        ax.set_thetagrids([0,90,180,270],labels=pos_label)
        ax.set_theta_direction("clockwise")
        ax.set_rorigin(-5)
        ax.set_theta_zero_location('N', offset=0)
    elif plot_type == "plane":
        ax = fig.add_subplot(111)
        for i, name_i in enumerate(name):
            c_i = color[i]
            r = cov_df[name_i]
            ax.plot(pos_ls, r, c=c_i, label=name_i)

        ax.set_ylabel(ylabel_str)
        ax.set_xlabel("mtDNA position")

    ax.legend()
    plt.tight_layout()
    fig.savefig(fig_name)
    plt.close()


def gen_mut_profile(bams, fmt="bam",sp="human",
             chr_name="MT",
             seq_type="sc",
             combined_bam=False,
             tag_name="CB",
             barcodes="NULL",
             min_read=200,
             max_depth=1e5,
             min_baseq=25,
             min_mapq=0,
             n_jobs=2):
    if isinstance(bams,str):
        sample = read_sam(bams,fmt,sp,chr_name,seq_type,combined_bam,
                          tag_name,barcodes,min_read,max_depth,min_baseq,min_mapq,n_jobs)
    elif isinstance(bams,Sample):
        sample=bams
    else:
        raise("ERROR: Invalid 'bams' value. "
              "It should be a path string to the folder of bam/sam files, or an result from function 'read_sam()'. ")

    mut_dict = {}
    mut_cov = {}
    index_ls_human = []
    index_ls_mus = []

    index_pos_human = []
    index_pos_mus = []
    base = ["A", "C", "G", "T"]
    for i in range(human_mt_len):
        for base_i in base:
            if human_mt_seq[i] != base_i:
                pos = i + 1
                index_ls_human.append(str(pos) + human_mt_seq[i] + ">" + base_i)
                index_pos_human.append(str(pos) + human_mt_seq[i] + ">" + base_i)

    for i in range(mus_mt_len):
        index_pos_mus.append(str(i + 1) + "_pos")
        for base_i in base:
            if mus_mt_seq[i] != base_i:
                pos = i + 1
                index_ls_mus.append(str(pos) + mus_mt_seq[i] + ">" + base_i)
                index_pos_mus.append(str(pos) + mus_mt_seq[i] + ">" + base_i)

    if isinstance(sample, list):
        for sample_i in sample:
            sp = sample_i.species
            if sample_i.seq_type == "bulk" or sample_i.seq_type == "sc":
                mut_dict[sample_i.name] = []
                mut_cov[sample_i.name] = []
                for i in range(sample_i.mt_len):
                    ref_i = sample_i.mtDNA[i]
                    if ref_i == "A":
                        AC_rate = float(sample_i.C[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        AG_rate = float(sample_i.G[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        AT_rate = float(sample_i.T[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        mut_dict[sample_i.name] += [AC_rate, AG_rate, AT_rate]
                    elif ref_i == "C":
                        CA_rate = float(sample_i.A[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        CG_rate = float(sample_i.G[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        CT_rate = float(sample_i.T[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        mut_dict[sample_i.name] += [CA_rate, CG_rate, CT_rate]
                    elif ref_i == "G":
                        GA_rate = float(sample_i.A[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        GC_rate = float(sample_i.C[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        GT_rate = float(sample_i.T[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        mut_dict[sample_i.name] += [GA_rate, GC_rate, GT_rate]
                    elif ref_i == "T":
                        TA_rate = float(sample_i.A[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        TC_rate = float(sample_i.C[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        TG_rate = float(sample_i.G[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        mut_dict[sample_i.name] += [TA_rate, TC_rate, TG_rate]
                    else:
                        NA_rate = float(sample_i.A[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        NC_rate = float(sample_i.C[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        NG_rate = float(sample_i.G[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        NT_rate = float(sample_i.T[i]) / sample_i.coverage[i] if sample_i.coverage[i] != 0 else 0.0
                        mut_dict[sample_i.name] += [NA_rate, NC_rate, NG_rate, NT_rate]

                    if ref_i == "N":
                        mut_cov[sample_i.name] += [sample_i.coverage[i], sample_i.coverage[i], sample_i.coverage[i],
                                                 sample_i.coverage[i]]
                    else:
                        mut_cov[sample_i.name] += [sample_i.coverage[i], sample_i.coverage[i], sample_i.coverage[i]]
            else:
                print("the data in sample "+sample_i.name+" is single cell sequencing data, please use the function sample.gen_mut_profile() for each sample.")
        if sp.lower() =="mus":
            profile_df = pd.DataFrame(data=mut_dict, index=index_ls_mus)
            df_coverage = pd.DataFrame(data=mut_cov, index=index_pos_mus)
        else:
            profile_df = pd.DataFrame(data=mut_dict, index=index_ls_human)
            df_coverage = pd.DataFrame(data=mut_cov, index=index_pos_human)

        adata = ad.AnnData(profile_df.T)
        adata.varm["coverage"] = df_coverage
        adata.uns["species"] = sp.lower()
        adata.uns['mut'] = True

        return adata
    elif isinstance(sample, Sample) and sample.seq_type=="sc":
        adata = sample.gen_mut_profile()
        return adata
    else:
        print("""Error (code 6): Wrong arguments for the function gen_mut_profile(). """)
        return -1

def add_metadata(adata,file_name,delimiter='\t'):
    df_metadata = pd.read_csv(file_name,sep=delimiter,index_col=0)
    if('label' not in df_metadata.columns):
        print("No column 'label' found in metadata, 'unknown' is used as the default cell labels")
        df_metadata['label'] = 'unknown'
    if('label_color' in df_metadata.columns):
        adata.uns['label_color'] = pd.Series(data=df_metadata.label_color.tolist(),index=df_metadata.label.tolist()).to_dict()
    else:
        print("No column 'label_color' found in metadata, random color is generated for each cell label")
        labels_unique = df_metadata['label'].unique()
        if(len(labels_unique)==1):
            adata.uns['label_color'] = {labels_unique[0]:'gray'}
        else:
            list_colors = sns.color_palette("hls",n_colors=len(labels_unique)).as_hex()
            adata.uns['label_color'] = {x:list_colors[i] for i,x in enumerate(labels_unique)}
        df_metadata['label_color'] = ''
        for x in labels_unique:
            id_cells = np.where(df_metadata['label']==x)[0]
            df_metadata.loc[df_metadata.index[id_cells],'label_color'] = adata.uns['label_color'][x]
    adata.obs = df_metadata.loc[adata.obs.index,:]
    adata.uns['metadata'] = True
    return None

def filter_features(adata, min_n_cells = None, max_n_cells=None,
                    min_pct_cells = None, max_pct_cells=None):

    if ('n_cells' in adata.var_keys()):
        n_cells = adata.var['n_cells']
    else:
        n_cells = np.sum(adata.X > 0, axis=0).astype(int)
        adata.var['n_cells'] = n_cells
    if ('pct_cells' in adata.var_keys()):
        pct_cells = adata.var['pct_cells']
    else:
        pct_cells = n_cells / adata.shape[0]
        adata.var['pct_cells'] = pct_cells

    if (sum(list(map(lambda x: x is None, [min_n_cells, min_pct_cells,
                                           max_n_cells, max_pct_cells]))) == 4):
        print('No filtering')
    else:
        feature_subset = np.ones(len(adata.var_names), dtype=bool)
        if (min_n_cells != None):
            print('Filtering based on min_n_cells')
            feature_subset = (n_cells >= min_n_cells) & feature_subset
        if (max_n_cells != None):
            print('Filtering based on max_n_cells')
            feature_subset = (n_cells <= max_n_cells) & feature_subset
        if (min_pct_cells != None):
            print('Filtering based on min_pct_cells')
            feature_subset = (pct_cells >= min_pct_cells) & feature_subset
        if (max_pct_cells != None):
            print('Filtering  based on max_pct_cells')
            feature_subset = (pct_cells <= max_pct_cells) & feature_subset
        adata._inplace_subset_var(feature_subset)
        adata.varm["coverage"] = adata.varm["coverage"].loc[adata.var_names,:]

        print("After filtering: "+str(adata.shape[0]) + ' cells, ' + str(adata.shape[1]) + ' features')
    return None


def filter_cells(adata,min_n_features=None, max_n_features=None,
                 min_pct_features=None, max_pct_features=None):

    if ('n_muts' in adata.obs_keys()):
        n_features = adata.obs['n_muts']
    else:
        n_features = np.sum(adata.X > 0, axis=1).astype(int)
        adata.obs['n_muts'] = n_features
    if ('pct_muts' in adata.obs_keys()):
        pct_features = adata.obs['pct_muts']
    else:
        pct_features = n_features / adata.shape[1]
        adata.obs['pct_muts'] = pct_features

    if (sum(list(map(lambda x: x is None, [min_n_features, min_pct_features,
                                           max_n_features, max_pct_features, ]))) == 4):
        print('No filtering')
    else:
        cell_subset = np.ones(len(adata.obs_names), dtype=bool)
        if (min_n_features != None):
            print('filter cells based on min_n_features')
            cell_subset = (n_features >= min_n_features) & cell_subset
        if (max_n_features != None):
            print('filter cells based on max_n_features')
            cell_subset = (n_features <= max_n_features) & cell_subset
        if (min_pct_features != None):
            print('filter cells based on min_pct_features')
            cell_subset = (pct_features >= min_pct_features) & cell_subset
        if (max_pct_features != None):
            print('filter cells based on max_pct_features')
            cell_subset = (pct_features <= max_pct_features) & cell_subset
        adata._inplace_subset_obs(cell_subset)

        adata.varm["coverage"] = adata.varm["coverage"].loc[:,adata.obs_names]

        print('after filtering : '+ str(adata.shape[0]) + ' cells, ' + str(adata.shape[1]) + 'mutations')
    return None

def plot_mean_std(adata, fig_name = "mean_std.png"):
    x_df = adata.to_df().T
    res_df = mut_mean_std(x_df)
    sns.set_theme(style="whitegrid")
    f, ax = plt.subplots(figsize=(6.5, 6.5))
    sns.despine(f, left=False, bottom=False)
    sns.scatterplot(x="mean", y="std",
                    palette="ch:r=-.2,d=.3_r",
                    s=10, linewidth=0,
                    data=res_df, ax=ax)
    # plt.title('Example Plot')
    # Set x-axis label
    plt.xlabel('Mean')
    # Set y-axis label
    plt.ylabel('Standard deviation')

    plt.tight_layout()
    f.savefig(fig_name)
    plt.close()

def cal_distance(adata,features="all"):
    if features=="all":
        af_df = adata.to_df()
        af_df_array = af_df.to_numpy()

        cov_df = adata.varm['coverage'].T
        cov_df_array = cov_df.to_numpy()
    elif features=="var_muts":
        af_df_array = adata.obsm['var_muts']

        cov_df = adata.varm['coverage'].T
        cov_df = cov_df.loc[adata.obs_names,adata.uns["var_muts"]]
        cov_df_array = cov_df.to_numpy()

    dis_df = pd.DataFrame(index=adata.obs_names,columns=adata.obs_names)
    n_sample = len(adata.obs_names)

    for index_i in range(n_sample):
        print(index_i)
        index_j = index_i
        while index_j < n_sample:
            af_ij = af_df_array[[index_i,index_j],:]
            cov_ij = cov_df_array[[index_i, index_j], :]

            Dij = cal_dis_ij(af_ij, cov_ij, c=0)

            dis_df.iloc[index_i,index_j] = Dij
            dis_df.iloc[index_j, index_i] = Dij
            index_j += 1
    adata.obsm["distance"] = dis_df
    return None


def select_specific_variants(adata,sample_name, min_af = 0.01, percent=0.8, check_other=True,other_min_af=0.01,other_n=5, other_percent=0):
    if isinstance(adata,ad.AnnData):
        x_df = adata.to_df().T
    else:
        x_df = adata
    x_df_sample = x_df.loc[:,x_df.columns.str.contains(sample_name)]

    x_df_selected = x_df_sample.loc[(x_df_sample > min_af).sum(axis=1)/len(x_df_sample.columns) > percent]
    var_ls = list(x_df_selected.index)

    if check_other:
        other_df = x_df.loc[:, ~x_df.columns.str.contains(sample_name)]
        other_df = other_df.loc[var_ls, :]
        percent_n = int(len(other_df.columns) * other_percent)
        other_n = other_n if percent_n < other_n else percent_n
        other_df_selected = other_df[(other_df > other_min_af).sum(axis=1) <= other_n]

        var_ls = list(other_df_selected.index)
    return var_ls


def select_top_varible_variants(adata, top_n = 50, percentile=20, fig_name = 'std_vs_means.png',fig_size = (8, 8),pad = 1.08):
    mean_muts = np.mean(adata.X, axis=0)
    std_muts = np.std(adata.X, ddof=1, axis=0)

    if (top_n is None):
        cutoff = np.percentile(std_muts, percentile)
        id_var_muts = np.where(std_muts > cutoff)[0]
        id_non_var_muts = np.where(std_muts <= cutoff)[0]
    else:
        id_var_muts = np.argsort(std_muts)[::-1][:top_n]
        id_non_var_muts = np.argsort(std_muts)[::-1][top_n:]

    adata.obsm['var_muts'] = adata.X[:, id_var_muts].copy()
    adata.uns['var_muts'] = adata.var_names[id_var_muts]


    ###plotting
    fig = plt.figure(figsize=fig_size)
    plt.scatter(mean_muts[id_var_muts], std_muts[id_var_muts], s=5, alpha=1, zorder=1, c='#EC4E4E')
    plt.scatter(mean_muts[id_non_var_muts], std_muts[id_non_var_muts], s=5, alpha=0.8, zorder=2, c='#6baed6')

    plt.xlabel('mean value')
    plt.ylabel('standard deviation')
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.tight_layout(pad=pad)

    plt.savefig(fig_name, pad_inches=1, bbox_inches='tight')
    plt.close(fig)
    return None


def select_top_principal_components(adata, features="var_muts", n_pc=50, max_pc=100, first_pc=True,fig_name='top_pcs.png',fig_size=(8, 8),
                                    pad=1.08, w_pad=None, h_pad=None):
    if (features is None):
        features = 'all'
    if (features not in ['all', 'var_muts']):
        raise ValueError("unrecognized feature '%s'" % features)

    sklearn_pca = PCA(svd_solver='arpack', random_state=12)
    if (features == 'var_muts'):
        print('using top variable muts ...')
        trans = sklearn_pca.fit(adata.obsm['var_muts'])
        X_pca = trans.transform(adata.obsm['var_muts'])
        pca_variance_ratio = trans.explained_variance_ratio_
        adata.obsm['pca'] = X_pca
        adata.uns['pca_variance_ratio'] = pca_variance_ratio
    else:
        print('using all the features ...')
        trans = sklearn_pca.fit(adata.X)
        X_pca = trans.transform(adata.X)
        pca_variance_ratio = trans.explained_variance_ratio_
        adata.obsm['pca'] = X_pca
        adata.uns['pca_variance_ratio'] = pca_variance_ratio

    adata.uns['top_pcs'] = trans

    max_pc = min(max_pc,len(pca_variance_ratio))

    if (first_pc):
        adata.obsm['top_pcs'] = X_pca[:, 0:(n_pc)]
    else:
        # discard the first Principal Component
        adata.obsm['top_pcs'] = X_pca[:, 1:(n_pc + 1)]
    print(str(n_pc) + ' PCs are selected')

    if ('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['select_top_principal_components'] = {'feature': features, 'first_pc': first_pc, 'n_pc': n_pc}

    ##plotting
    fig = plt.figure(figsize=fig_size)
    plt.plot(range(max_pc), pca_variance_ratio[:max_pc])
    if (first_pc):
        plt.axvline(n_pc, c='red', ls='--')
    else:
        plt.axvline(1, c='red', ls='--')
        plt.axvline(n_pc + 1, c='red', ls='--')
    plt.xlabel('Principal Component')
    plt.ylabel('Variance Ratio')
    plt.locator_params(axis='x', nbins=5)
    plt.locator_params(axis='y', nbins=5)
    plt.tight_layout(pad=pad, h_pad=h_pad, w_pad=w_pad)

    plt.savefig(fig_name, pad_inches=1, bbox_inches='tight')
    plt.close(fig)
    return None

def plot_clustermap(adata, x_df=None, features="var_muts", ann_color=None, ann_label=None, fig_name="clusterheatmap.png", vmin=0, vmax=1,
                    fig_size=(12,12), xlabel=False,ylabel=True,xlabel_size=8,ylabel_size=8,xcluster=True, ycluster=True):
    if not('metadata' in adata.uns and adata.uns['metadata']):
        ann_color=None

    if isinstance(x_df,pd.DataFrame):
        x_df = x_df
    else:
        if features=="all":
            x_df = adata.to_df().T
            x_df.fillna(value=np.nan, inplace=True)
        elif features=="var_muts":
            if "var_muts" in adata.obsm_keys():
                x_array = adata.obsm['var_muts']
                x_df = pd.DataFrame(data=x_array, index = adata.obs_names,columns=adata.uns["var_muts"])
                x_df = x_df.T
                x_df.fillna(value=np.nan, inplace=True)
            else:
                print("Fearures '"+features+"' were selected to use. Runn ‘select_top_varible_variants’ first to set top N variable mutation positions.")
                return -1
        elif features=="distance":
            if "distance" in adata.obsm_keys():
                x_df = adata.obsm['distance']
                x_df.fillna(value=np.nan, inplace=True)
            else:
                print("Running ‘cal_distance’ first to calculate the cell distance matrix.")
                return -1
        else:
            print("Invalid value for 'features', it should be in ['all','var_muts'(run select_top_varible_variants() first),'distance'(run cal_distance() first)].")
            return -1

    obs_df = adata.obs.loc[x_df.columns,:]

    color_pal_ls = ["hls", "tab10", "Spectral", "coolwarm", "Paired", "YlOrBr"]

    color2label={}
    ann_color_dict = {}
    if ann_color == None:
        col_colors_df=None
    elif isinstance(ann_color,list):
        if isinstance(ann_label, list) and len(ann_color) == len(ann_label):
            pass
        else:
            ann_label = ann_color
            print("WARNING: 'ann_label' should be a list with same length as 'ann_color'. 'ann_label' was ignored!'")
        for i, ann_i in enumerate(ann_color):
            if ann_i in obs_df.columns:
                labels_unique = obs_df[ann_i].unique()
                if (len(labels_unique) == 1):
                    adata.uns[ann_i+"_color"] = {labels_unique[0]: 'gray'}
                else:
                    list_colors = sns.color_palette(color_pal_ls[i%6], n_colors=len(labels_unique)).as_hex()
                    adata.uns[ann_i+"_color"] = {x: list_colors[i] for i, x in enumerate(labels_unique)}
                col_colors = list(pd.Series(obs_df[ann_i]).map(adata.uns[ann_i+"_color"]))
                col_colors_name = ann_label[i]
                ann_color_dict[col_colors_name] = col_colors
                color2label[ann_i] = col_colors_name
            else:
                print("WARNING: There is no '" + ann_i + "' annotation in the anndata.obs_names. Ignored!'")

        col_colors_df = pd.DataFrame(data=ann_color_dict,index=x_df.columns)
    else:
        print("WARNING: 'ann_color' should be a list witn elements in the metadata file")

    c = [ "white", "lightcoral", "red", "darkred"]
    v = [0, 0.33, 0.67, 1.]
    l = list(zip(v, c))
    cmap = LinearSegmentedColormap.from_list('wr', l, N=256)
    sns.set_theme()

    g=sns.clustermap(x_df, cmap=cmap,vmin=vmin, vmax=vmax,xticklabels=xlabel,yticklabels=ylabel,
                     row_cluster=xcluster, col_cluster=xcluster,
                     figsize=fig_size, col_colors = col_colors_df)
    g.cax.set_position([.97, .25, .03, .45])

    # The following two lines generate custom fake lines that will be used as legend entries:
    pos_i = 0
    l={}
    for color_i in color2label:
        markers = [plt.Line2D([0, 0], [0, 0], color=color, marker='o', linestyle='') for color in adata.uns[color_i+"_color"].values()]
        l[pos_i]=plt.legend(markers, adata.uns[color_i+"_color"].keys(), title= color2label[color_i], ncol=1,
                   numpoints=1,loc="upper left", bbox_to_anchor=(2.1, 1.0 - pos_i), frameon=False)
        pos_i +=  0.15
        pos_i += 0.035 * (len(adata.uns[color_i + "_color"].keys())+1)
    for pos_i in l:
        plt.gca().add_artist(l[pos_i])

    if xlabel:
        g.ax_heatmap.tick_params(axis='x', labelsize=xlabel_size)
    if ylabel:
        g.ax_heatmap.tick_params(axis='y', labelsize=ylabel_size)

    if not xcluster:
        g.ax_col_dendrogram.remove()
    if not ycluster:
        g.ax_row_dendrogram.remove()

    plt.tight_layout()
    g.savefig(fig_name)
    plt.close()


def dimension_reduction(adata,features="all", n_neighbors=20, n_components=2, dim_reducer="tsne", **kwargs):
    if features == "all":
        x_df = adata.to_df()
    elif features == "var_muts":
        if "var_muts" in adata.obsm_keys():
            x_array = adata.obsm['var_muts']
            x_df = pd.DataFrame(data=x_array, index=adata.obs_names, columns=adata.uns["var_muts"])
            x_df.fillna(value=np.nan, inplace=True)
        else:
            print("Running ‘select_variant’ first to set top N variable mutation positions.")
            return -1
    elif features == "distance":
        if "distance" in adata.obsm_keys():
            x_df = adata.obsm['distance']
            x_df.fillna(value=np.nan, inplace=True)
        else:
            print("Running ‘cal_distance’ first to calculate the cell distance matrix.")
            return -1
    elif features == "top_pcs":
        if "top_pcs" in adata.obsm_keys():
            x_df = adata.obsm['top_pcs']
        else:
            print("Running ‘select_top_principal_components’ first to calculate the cell distance matrix.")
            return -1

    else:
        print(
            "Invalid value for 'feature', it should be in ['all','var_muts'(run select_variant() first),'distance'(run cal_distance() first)].")
        return -1

    dim_reducer=dim_reducer.upper()
    if dim_reducer == "TSNE":
        tsne = TSNE(n_components=n_components, perplexity=30.0, n_iter=1000,random_state=12)
        trans = tsne.fit(x_df)
        adata.uns['trans_tsne'] = trans
        adata.obsm['X_tsen'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    elif dim_reducer == "UMAP":
        reducer = umap.UMAP(n_neighbors=n_neighbors, n_components=n_components, random_state=12)
        trans = reducer.fit(x_df)
        adata.uns['trans_umap'] = trans
        adata.obsm['X_umap'] = trans.embedding_
        adata.obsm['X_dr'] = trans.embedding_
    elif dim_reducer == "PCA":
        reducer = PCA(n_components=n_components,svd_solver='arpack',random_state=12)
        trans = reducer.fit(x_df)
        adata.uns['trans_pca'] = trans
        adata.obsm['X_pca'] = trans.transform(x_df)
        adata.obsm['X_dr'] = adata.obsm['X_pca']
    else:
        print("Error (code 7): Wrong value for argument 'method' ('tSNE','UMAP','PCA').")
        return -1

    if ('params' not in adata.uns_keys()):
        adata.uns['params'] = dict()
    adata.uns['params']['dimension_reduction'] = {'feature': features, 'method': dim_reducer, 'n_components': n_components}
    return None


def plot_clusters(adata, fig_name="cluster.png", cluster_color="auto",method="DBSCAN",cluster_n=2,color_pal="Blues",cbar_frac=0.025, **kwargs):
    if 'params' not in adata.uns_keys() or 'dimension_reduction' not in adata.uns['params']:
        print("ERROR: Please run 'dimension_reduction' first.")
        return -1
    dimension =  adata.uns['params']['dimension_reduction']['n_components']
    dim_reducer = adata.uns['params']['dimension_reduction']['method']

    X_embedded = adata.obsm['X_dr']

    if cluster_color == "auto":
        print("Cluster color is set to 'auto'. Sample cluster will be detected using algorithm " + method)
        if method.upper() == "KMEANS":
            cluster_labels = KMeans(n_clusters=cluster_n, random_state=0).fit_predict(X_embedded)
        else:
            cluster_labels = DBSCAN(eps=3, min_samples=2).fit_predict(X_embedded)
        adata.obs["detect_cluster"] = cluster_labels
        ylabel = np.asarray(cluster_labels)
    else:
        if cluster_color not in adata.obs.columns:
            print("Invalid value of 'cluster_color', it should be one of the values in anndata.obs_names.")
            return -1
        else:
            ylabel = adata.obs[cluster_color]

    y_set = list(set(ylabel))
    color_len = len(y_set)
    palette_ls = sns.husl_palette(color_len)
    if color_len < 32:
        if dimension == 2:
            sns.set(style="ticks")
            f, ax = plt.subplots(figsize=(12, 9))
            sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=ylabel,legend="full", palette=palette_ls,ax =ax)
            dim_reducer = dim_reducer if dim_reducer != "TSNE" else "tSNE"
            ax.set(xlabel=dim_reducer + " 1")
            ax.set(ylabel=dim_reducer + " 2")
            ax.legend(loc=2, bbox_to_anchor=(1, 1), frameon=False)

            plt.tight_layout()
            f.savefig(fig_name)
            plt.close()

        if dimension == 3:
            fig = plt.figure(figsize=(12, 9))
            ax = fig.add_subplot(111, projection='3d')
            for y_i in range(len(y_set)):
                ax.scatter(X_embedded[np.where(ylabel == y_set[y_i]), 0],
                           X_embedded[np.where(ylabel == y_set[y_i]), 1],
                           X_embedded[np.where(ylabel == y_set[y_i]), 2],
                           s=10, c=palette_ls[y_i], label=y_set[y_i])

            dim_reducer = dim_reducer if dim_reducer != "TSNE" else "tSNE"
            ax.set_xlabel(dim_reducer + " 1")
            ax.set_ylabel(dim_reducer + " 2")
            ax.set_zlabel(dim_reducer + " 3")
            ax.legend()
            fig.savefig(fig_name)
            plt.close()
    else:
        palette_ls = sns.color_palette(color_pal, as_cmap=True)

        norm = plt.Normalize(ylabel.min(), ylabel.max())
        sm = plt.cm.ScalarMappable(cmap=palette_ls, norm=norm)
        sm.set_array([])

        if dimension == 2:
            sns.set(style="ticks")
            f, ax = plt.subplots(figsize=(12, 9))
            g = sns.scatterplot(x=X_embedded[:, 0], y=X_embedded[:, 1], hue=ylabel, palette=palette_ls,ax=ax)

            m0 = int(np.floor(ylabel.min()))  # colorbar min value
            m4 = int(np.ceil(ylabel.max()))  # colorbar max value
            m1 = int(1 * (m4 - m0) / 4.0 + m0)  # colorbar mid value 1
            m2 = int(2 * (m4 - m0) / 4.0 + m0)  # colorbar mid value 2
            m3 = int(3 * (m4 - m0) / 4.0 + m0)  # colorbar mid value 3
            # use this colorbar
            cbar = plt.colorbar(sm, ax=g.axes)
            cbar.set_ticks([m0, m1, m2, m3, m4])
            cbar.set_ticklabels([m0, m1, m2, m3, m4])

            dim_reducer = dim_reducer if dim_reducer != "TSNE" else "tSNE"
            ax.set(xlabel=dim_reducer + " 1")
            ax.set(ylabel=dim_reducer + " 2")
            ax.legend(loc=2,bbox_to_anchor=(1.15, 1),frameon=False)
            plt.tight_layout()
            f.savefig(fig_name)
            plt.close()

        if dimension == 3:
            fig = plt.figure(figsize=(12, 9))
            ax = fig.add_subplot(111, projection='3d')

            ax.scatter(X_embedded[:, 0],X_embedded[:, 1],X_embedded[:,2],
                           s=10, cmap=palette_ls, c=ylabel)

            m0 = int(np.floor(ylabel.min()))  # colorbar min value
            m4 = int(np.ceil(ylabel.max()))  # colorbar max value
            m1 = int(1 * (m4 - m0) / 4.0 + m0)  # colorbar mid value 1
            m2 = int(2 * (m4 - m0) / 4.0 + m0)  # colorbar mid value 2
            m3 = int(3 * (m4 - m0) / 4.0 + m0)  # colorbar mid value 3
            # use this colorbar
            cbar = plt.colorbar(sm, ax=ax,fraction=cbar_frac, pad=0.04)
            cbar.set_ticks([m0, m1, m2, m3, m4])
            cbar.set_ticklabels([m0, m1, m2, m3, m4])

            dim_reducer = dim_reducer if dim_reducer != "TSNE" else "tSNE"
            ax.set_xlabel(dim_reducer + " 1")
            ax.set_ylabel(dim_reducer + " 2")
            ax.set_zlabel(dim_reducer + " 3")
            # ax.legend()

            fig.savefig(fig_name)
            plt.close()

