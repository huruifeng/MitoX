from .global_data import human_mt_seq
from .global_data import human_mt_len

from .global_data import mus_mt_seq
from .global_data import mus_mt_len
from .bamreader import *

import os
import shutil

import pandas as pd
import numpy as np
import anndata as ad

from multiprocessing import Pool
from datetime import datetime

import matplotlib.pyplot as plt
import seaborn as sns

import HTSeq
import pysam

pysam.set_verbosity(0)

import tempfile

# cmd_str = "python -m HTSeq.scripts.count -i ENSG00000198727 --additional-attr=gene_name -f bam bam/AZ_A2.bam Homo_sapiens.GRCh38.101.gtf.gz >> out.txt"

human_gene_id = { "MT-ATP8": "ENSG00000228253","MT-ATP6": "ENSG00000198899","MT-CO1": "ENSG00000198804","MT-CO2": "ENSG00000198712",
            "MT-CO3": "ENSG00000198938", "MT-CYB": "ENSG00000198727", "MT-ND1": "ENSG00000198888","MT-ND2": "ENSG00000198763",
            "MT-ND3": "ENSG00000198840", "MT-ND4L": "ENSG00000212907","MT-ND4": "ENSG00000198886", "MT-ND5": "ENSG00000198786",
            "MT-ND6": "ENSG00000198695", "MT-TA": "ENSG00000210127","MT-TR": "ENSG00000210174", "MT-TN": "ENSG00000210135",
            "MT-TD": "ENSG00000210154",  "MT-TC": "ENSG00000210140","MT-TE": "ENSG00000210194","MT-TQ": "ENSG00000210107",
            "MT-TG": "ENSG00000210164", "MT-TH": "ENSG00000210176","MT-TI": "ENSG00000210100","MT-TL1": "ENSG00000209082",
            "MT-TL2": "ENSG00000210191","MT-TK": "ENSG00000210156","MT-TM": "ENSG00000210112","MT-TF": "ENSG00000210049",
            "MT-TP": "ENSG00000210196", "MT-TS1": "ENSG00000210151","MT-TS2": "ENSG00000210184","MT-TT": "ENSG00000210195",
            "MT-TW": "ENSG00000210117","MT-TY": "ENSG00000210144","MT-TV": "ENSG00000210077","MT-RNR1": "ENSG00000211459",
            "MT-RNR2": "ENSG00000210082"}
human_id_gene = {v: k for k, v in human_gene_id.items()}

mus_gene_id = {"mt-Atp6":"ENSMUSG00000064357","mt-Atp8":"ENSMUSG00000064356","mt-Co1":"ENSMUSG00000064351","mt-Co2":"ENSMUSG00000064354",
              "mt-Co3":"ENSMUSG00000064358","mt-Cytb":"ENSMUSG00000064370","mt-Nd1":"ENSMUSG00000064341","mt-Nd2":"ENSMUSG00000064345",
              "mt-Nd3":"ENSMUSG00000064360","mt-Nd4":"ENSMUSG00000064363","mt-Nd4l":"ENSMUSG00000065947","mt-Nd5":"ENSMUSG00000064367",
              "mt-Nd6":"ENSMUSG00000064368","mt-Rnr1":"ENSMUSG00000064337","mt-Rnr2":"ENSMUSG00000064339","mt-Ta":"ENSMUSG00000064347",
              "mt-Tc":"ENSMUSG00000064349","mt-Td":"ENSMUSG00000064353","mt-Te":"ENSMUSG00000064369","mt-Tf":"ENSMUSG00000064336",
              "mt-Tg":"ENSMUSG00000064359","mt-Th":"ENSMUSG00000064364","mt-Ti":"ENSMUSG00000064342","mt-Tk":"ENSMUSG00000064355",
              "mt-Tl1":"ENSMUSG00000064340","mt-Tl2":"ENSMUSG00000064366","mt-Tm":"ENSMUSG00000064344","mt-Tn":"ENSMUSG00000064348",
              "mt-Tp":"ENSMUSG00000064372","mt-Tq":"ENSMUSG00000064343","mt-Tr":"ENSMUSG00000064361","mt-Ts1":"ENSMUSG00000064352",
              "mt-Ts2":"ENSMUSG00000064365","mt-Tt":"ENSMUSG00000064371","mt-Tv":"ENSMUSG00000064338","mt-Tw":"ENSMUSG00000064346",
              "mt-Ty":"ENSMUSG00000064350"}
mus_id_gene = {v: k for k, v in mus_gene_id.items()}

def do_reads_count(params):
    (bam_folder, samfile, i, sam_n, expression_temp, chr_name, seqlen, min_mapq, exons) = params

    now = datetime.now()  # current date and time
    date_time = now.strftime("%b/%d/%Y %a %H:%M:%S")
    print(date_time + " Processing " + samfile + "... [" + str(i) + "/" + str(sam_n) + "]\n")

    sam_aln = pysam.AlignmentFile(bam_folder + "/" + samfile)
    if not sam_aln.has_index():
        print("No index file was found. Indexing the BAM file...")
        pysam.index(bam_folder + "/" + samfile)
    sam_aln.close()

    mt_sam_name = "mt_" + samfile[:-4]
    bamfile = pysam.AlignmentFile(bam_folder + "/" + samfile)
    mtreads = pysam.AlignmentFile(expression_temp + "/" + mt_sam_name + ".bam", "wb", template=bamfile)
    for read in bamfile.fetch(contig=chr_name, start=0, stop=seqlen):
        if read.mapping_quality > min_mapq:
            mtreads.write(read)
    bamfile.close()
    mtreads.close()
    pysam.index(expression_temp + "/" + mt_sam_name + ".bam")

    counts = {}
    bam_file = HTSeq.BAM_Reader(expression_temp + "/" + mt_sam_name + ".bam")
    for alnmt in bam_file:
        if alnmt.aligned:
            iset = None
            for iv2, step_set in exons[alnmt.iv].steps():
                if iset is None:
                    iset = step_set.copy()
                else:
                    iset.intersection_update(step_set)
            if len(iset) == 1:
                if list(iset)[0] in counts:
                    counts[list(iset)[0]] += 1
                else:
                    counts[list(iset)[0]] = 1
    return {samfile[:-4]: counts}

def gen_expression_profile(bam_folder, gtf_file, fmt="bam", sp="human",
                           chr_name="MT",
                           seq_type="sc",
                           combined_bam=False,
                            tag_name="CB",
                            barcodes="NULL",
                            min_read=200,
                           stranded=True,
                           feature_type="exon",
                           min_mapq=10,
                           n_jobs=2):

    folder = bam_folder
    files = os.listdir(folder)
    sam_ls = []

    if fmt == "bam" or fmt == "sam":
        if fmt == "bam":
            rb = "rb"
            for file_i in files:
                if file_i.endswith(".bam") and os.path.isfile(folder + "/" + file_i):
                    sam_ls.append(file_i)
        else:
            rb = "r"
            for file_i in files:
                if file_i.endswith(".sam") and os.path.isfile(folder + "/" + file_i):
                    sam_ls.append(file_i)
    else:
        raise Exception("""Error (Code 1): Parameter 'fmt' setting is wrong (Valid value: 'bam' OR 'sam').""")
        return -1
    sam_n = len(sam_ls)

    seq_len = mus_mt_len if sp.lower() == "mus" else human_mt_len

    ######
    if seq_type == "bulk":
        read_adata = sepearted_bam_profile(bam_folder, gtf_file, fmt, sp, chr_name, stranded,feature_type, min_mapq,n_jobs)
    elif seq_type == "sc":
        if not combined_bam:
            print("Each BAM file are considered as a cell sample and running on each BAM file separately...")
            read_adata = sepearted_bam_profile(bam_folder, gtf_file, fmt, sp, chr_name, stranded, feature_type,
                                               min_mapq, n_jobs)
        else:
            print("Each BAM file contains many cells, and barcodes will be used for cell count...")
            print(str(n_jobs) + ' processes are running...')

            now = datetime.now()  # current date and time
            date_time = now.strftime("%b/%d/%Y %a %H:%M:%S")
            print(date_time + " Reading BAM/SAM files...")

            ##
            params = [(folder, sam_ls[i], rb) for i in range(sam_n)]
            pool = Pool(processes=n_jobs)
            results = pool.map(do_index, params)
            pool.close()

            ##
            good_barcodes = {}
            barcodes_in_bam = {}
            params = [(folder, sam_ls[i], rb, chr_name, 0, seq_len, tag_name, barcodes, min_read) for i in range(sam_n)]
            pool = Pool(processes=n_jobs)
            results = pool.map(do_get_barcodes, params)
            pool.close()
            for x_i in results:
                good_barcodes[x_i[0]] = x_i[1]
                barcodes_in_bam[x_i[0]] = x_i[2]

            ##
            tmp_split_folder = {}
            params = [(folder, sam_ls[i], rb, good_barcodes[sam_ls[i]], tag_name, chr_name, n_jobs) for i in
                      range(sam_n)]
            pool = Pool(processes=n_jobs)
            results = pool.map(do_split_bam, params)
            pool.close()

            for x_i in results:
                tmp_split_folder[x_i[0]] = x_i[1]

            ######
            expr_ls = []
            for sam_i in sam_ls:
                fmt="bam"
                print("Processing file: "+sam_i)
                expr_adata = sepearted_bam_profile(tmp_split_folder[sam_i], gtf_file, fmt, sp, chr_name, stranded, feature_type,
                                                   min_mapq, n_jobs)

                expr_ls.append(expr_adata.to_df().T)

            expr_df = pd.concat(expr_ls,axis=1)
            read_adata = ad.AnnData(expr_df)
    return read_adata


def sepearted_bam_profile(bam_folder, gtf_file, fmt="bam", sp="human", chr_name="MT", stranded=True,
                           feature_type="exon", min_mapq=10, n_jobs=2):

    if sp.lower() == "human":
        seqlen = human_mt_len
    elif sp.lower() == "mus":
        seqlen = mus_mt_len
    else:
        raise ("Error: Invalid 'sp' value. It should be 'human' or 'mus'.")

    folder = bam_folder
    files = os.listdir(folder)
    sam_ls = []
    if fmt not in ["bam", "sam"]:
        raise Exception("""Error (Code 11): Parameter 'fmt' setting is wrong (Valid value: 'bam' OR 'sam').""")

    for file_i in files:
        if file_i.endswith(".bam") and os.path.isfile(folder + "/" + file_i):
            sam_ls.append(file_i)
    sam_n = len(sam_ls)
    print(str(sam_n) + " BAM files were detected.")

    print("Preparing GTF file...")
    gtf_anno = HTSeq.GFF_Reader(gtf_file)
    exons = HTSeq.GenomicArrayOfSets("auto", stranded=stranded)

    for feature in gtf_anno:
        if feature.iv.chrom != "MT":
            continue
        if feature.type == feature_type:
            exons[feature.iv] += feature.name

    expression_temp = tempfile.mkdtemp()
    params = [(bam_folder, sam_ls[i], i, sam_n, expression_temp, chr_name, seqlen, min_mapq, exons) for i in
              range(sam_n)]
    ##
    pool = Pool(processes=n_jobs)
    results = pool.map(do_reads_count, params)
    pool.close()

    if sp.lower()=="human":
        gene_id = human_gene_id
        id_gene = human_id_gene
        start_str = "ENSG"
    elif sp.lower()=="mus":
        gene_id = mus_gene_id
        id_gene = mus_id_gene
        start_str = "ENSMU"
    else:
        raise ("ERROR: Invalid 'sp' value, it should be 'human' or 'mus'.")

    read_counts = {}
    for i in range(sam_n):
        sample_name = list(results[i].keys())[0]
        read_counts[sample_name] = {}
        gene_ls = results[i][sample_name]
        id_x = list(gene_ls.keys())[0]

        if id_x.startswith(start_str):
            for gene_i in id_gene:
                if gene_i in gene_ls.keys():
                    read_counts[sample_name][id_gene[gene_i]] = gene_ls[gene_i]
                else:
                    read_counts[sample_name][id_gene[gene_i]] = 0
        else:
            for gene_i in gene_id:
                if gene_i in gene_ls.keys():
                    read_counts[sample_name][gene_i] = gene_ls[gene_i]
                else:
                    read_counts[sample_name][gene_i] = 0

    print("Clearning temporary files...")
    shutil.rmtree(expression_temp)

    read_df = pd.DataFrame(data=read_counts)
    read_adata = ad.AnnData(read_df.T)
    read_adata.uns["expr"] = True
    return read_adata


def combine_mut_expr(mut_adata, expr_adata):
    combined_adata = mut_adata.copy()
    df_expr = expr_adata.to_df()

    x_1 = combined_adata.obs_names
    x_2 = df_expr.index
    x = list(set(x_1) & set(x_2))

    df_expr = df_expr.loc[x, :]
    combined_adata = combined_adata[x, :]



    if "metadata" in mut_adata.uns_keys():
        pass
    elif "metadata" in expr_adata.uns_keys():
        for meta_i in expr_adata.obs_keys():
            combined_adata.obs[meta_i] = expr_adata.obs.loc[x,meta_i]
    else:
        raise("ERROR: No metadata was found, please add metadata to the mut_adata or expr_adata first..")

    for col_i in df_expr.columns:
        combined_adata.obs[col_i] = df_expr[col_i]

    combined_adata.uns['combine'] = True
    return combined_adata

def plot_expr(adata, gene, group=None, using={}, fig_name="expr_plot.png", fig_size=(12, 8),
              plot_type="boxplot",strip=False,log2=True, all_in_one=True,col_n=0,xrotate=30, **kwargs):
    if not(isinstance(gene,list) and len(gene)) > 0:
        raise ("Please set a list of genes.")

    if "expr" in adata.uns and adata.uns["expr"]:
        expr_df = adata.to_df()
        for gene_i in gene:
            if gene_i not in expr_df.columns:
                print("WARNING: "+gene_i+ " is not in the data, please check if the name is correct.")
                gene.remove(gene_i)
        ylabel="Expression (read counts)"
        if log2:
            expr_df = np.log2(expr_df+1)
            ylabel="Expression (log2(counts+1))"

        if not ("metadata" in adata.uns and adata.uns["metadata"]):
            expr_gene_df = expr_df.loc[:, gene]
            data_df = expr_gene_df.melt(var_name="gene")
            print("WARNING: There is no annotation data for samples.'group' was set to None .")
            group = None
        else:
            df_meta = adata.obs.copy(deep=True)
            meta_cols = list(df_meta.columns)
            expr_meta_df = pd.concat([expr_df, df_meta], axis=1, join='inner')
            gene_and_meta_col = gene + meta_cols
            expr_meta_df = expr_meta_df.loc[:, gene_and_meta_col]
            data_df = expr_meta_df.melt(id_vars=meta_cols, var_name="gene")

        if len(using) != 0:
            k_col = list(using.keys())[0]
            k_values = using[k_col]
            data_df = data_df.loc[data_df[k_col].isin(k_values),:]

    elif "combine" in adata.uns and adata.uns["combine"]:
        # expr_df = adata.to_df()
        for gene_i in gene:
            if gene_i not in adata.obs_names:
                print("WARNING: " + gene_i + " is not in the data, please check if the name is correct.")
                gene.remove(gene_i)
        ylabel = "Expression (read counts)"

        expr_df = adata.obs.loc[:,gene]

        if log2:
            expr_df = np.log2(expr_df + 1)
            ylabel = "Expression (log2(counts+1))"

        if not ("metadata" in adata.uns and adata.uns["metadata"]):
            expr_gene_df = expr_df.loc[:, gene]
            data_df = expr_gene_df.melt(var_name="gene")
            print("WARNING: There is no annotation data for samples.'group' was set to None .")
            group = None
        else:
            df_meta = adata.obs.copy(deep=True)
            for gene_i in gene:
                if gene_i in df_meta.columns:
                    df_meta.drop(gene_i,axis=1,inplace=True)

            meta_cols = list(df_meta.columns)
            expr_meta_df = pd.concat([expr_df, df_meta], axis=1, join='inner')
            gene_and_meta_col = gene + meta_cols
            expr_meta_df = expr_meta_df.loc[:, gene_and_meta_col]
            data_df = expr_meta_df.melt(id_vars=meta_cols, var_name="gene")

        if len(using) != 0:
            k_col = list(using.keys())[0]
            k_values = using[k_col]
            data_df = data_df.loc[data_df[k_col].isin(k_values), :]
    else:
        raise("ERROR: The 'adata' should should be expression anndata or combined anndata.")

    ###########
    fig = plt.figure(num=1, figsize=fig_size, dpi=80)
    if all_in_one:
        ax = fig.add_subplot(111)
        sns.set_style(style="ticks")
        if group == None:
            if plot_type == "boxplot":
                sns.boxplot(x="gene", y="value", data=data_df, dodge=True, **kwargs)
            elif plot_type == "violin":
                sns.violinplot(x="gene", y="value", data=data_df, dodge=True, **kwargs)
            if strip:
                sns.stripplot(x="gene", y="value", data=data_df, dodge=True, **kwargs)

        else:
            if plot_type == "boxplot":
                sns.boxplot(x="gene", y="value", hue=group, data=data_df, dodge=True, **kwargs)
            elif plot_type == "violin":
                sns.violinplot(x="gene", y="value", hue=group, data=data_df, dodge=True, **kwargs)

            if strip:
                sns.stripplot(x="gene", y="value", hue=group, data=data_df, dodge=True, **kwargs)

            n_color = len(data_df[group].unique())
            handles, labels = ax.get_legend_handles_labels()
            l = plt.legend(handles[0:n_color], labels[0:n_color], title=group, bbox_to_anchor=(1.05, 1), loc=2,
                           borderaxespad=0., frameon=False)
        ax.set_title("")
        ax.set_xlabel("")
        ax.set_ylabel(ylabel)
        sns.despine(offset=10, top=True, right=True)
        fig.tight_layout()
    else:
        gene_n = len(gene)
        col_n = gene_n if col_n == 0 else col_n
        if gene_n % col_n == 0:
            row_n = int(gene_n / col_n)
        else:
            row_n = int(gene_n / col_n) + 1
        i = 0
        ii = 0
        ax_ls = []
        for i in range(gene_n):
            i += 1
            gene_ii = gene[ii]
            ii += 1
            print(gene_ii)
            data_gene_df = data_df.loc[data_df["gene"] == gene_ii, :]
            ax = fig.add_subplot(row_n, col_n, i)
            ax_ls.append(ax)
            # print(row_n,col_n,i)
            sns.set_style(style="ticks")
            if group == None:
                if plot_type == "boxplot":
                    g = sns.boxplot(x="gene", y="value", data=data_gene_df, dodge=True, **kwargs)
                elif plot_type == "violin":
                    g = sns.violinplot(x="gene", y="value", data=data_gene_df, dodge=True, **kwargs)
                if strip:
                    g = sns.stripplot(x="gene", y="value", data=data_gene_df, dodge=True, **kwargs)

            else:
                if plot_type == "boxplot":
                    g = sns.boxplot(x=group, y="value", data=data_gene_df, **kwargs)
                elif plot_type == "violin":
                    g = sns.violinplot(x=group, y="value", data=data_gene_df, **kwargs)

                if strip:
                    g = sns.stripplot(x=group, y="value", data=data_gene_df, **kwargs)
            ax.set_title(gene_ii)
            ax.set_xlabel("")
            ax.set_ylabel("")
            sns.despine(offset=10, top=True, right=True)
            i = i % col_n

        ax_all = []
        for ax_i in ax_ls:
            ax_all += ax_i.get_xticklabels()
        plt.setp(ax_all, rotation=xrotate)

    fig.tight_layout()
    fig.savefig(fig_name)
    plt.close(fig)



