import os
import shutil

import pysam
pysam.set_verbosity(0)
from datetime import datetime
from multiprocessing import Pool

import tempfile


from .sample import Sample, Cell
from .bamwriter import split_bam

from .global_data import human_mt_seq
from .global_data import human_mt_len

from .global_data import mus_mt_seq
from .global_data import mus_mt_len


def test():
    print("Test done.")


def do_coverage_count(params):
    # print(params)
    (folder, samfile,i,sam_n, rb, chr_name, start_x, stop_x, min_baseq_x, min_mapq_x) = params

    ## read_callback function
    def keep_read(read):
        if (read.flag & (0x4 | 0x100 | 0x200 | 0x400)) or read.mapping_quality < min_mapq_x:
            return False
        else:
            return True

    now = datetime.now()  # current date and time
    date_time = now.strftime("%b/%d/%Y %a %H:%M:%S")
    print(date_time + " Processing " + samfile + "... ["+str(i)+"/"+str(sam_n)+"]\n")

    bam_path = folder + "/" + samfile
    sam_aln = pysam.AlignmentFile(bam_path, rb)
    if not sam_aln.has_index():
        print("No index file was found. Indexing the BAM file...")
        pysam.index(bam_path)
    sam_aln.close()

    sam_aln = pysam.AlignmentFile(bam_path, rb)
    ACGT_cov_x = sam_aln.count_coverage(chr_name, start=start_x, stop=stop_x, quality_threshold=min_baseq_x,
                                        read_callback=keep_read)

    reads_n = sam_aln.count(chr_name, start=start_x, stop=stop_x)

    # MTpileup = samfile.pileup('MT', start=0,stop=seq_len, stepper='all', max_depth=500000, truncate=True,
    #                           min_base_quality=min_baseq, min_mapping_quality=min_mapq,ignore_overlaps=False)
    sam_aln.close()

    coverage = []
    (A,C,G,T) = (ACGT_cov_x[0], ACGT_cov_x[1], ACGT_cov_x[2], ACGT_cov_x[3])
    for i in range(start_x,stop_x):
        coverage_i = A[i] + C[i] + G[i] + T[i]
        coverage.append(coverage_i)

    now = datetime.now()  # current date and time
    date_time = now.strftime("%b/%d/%Y %a %H:%M:%S")
    print(date_time + " Processing " + samfile + "... Done !\n")
    return (A,C,G,T,coverage,reads_n)

def do_index(params):
    (folder, samfile, rb) = params
    bam_path = folder + "/" + samfile
    sam_aln = pysam.AlignmentFile(bam_path, rb)
    if not sam_aln.has_index():
        print("["+samfile+"] No index file was found. Indexing the BAM file...")
        pysam.index(bam_path)
        print("[" + samfile + "] Indexing complete...")
    sam_aln.close()
    return 1

def do_get_barcodes(params):
    (folder, samfile, rb, chr_name, start_x,stop_x ,tag_name, barcodes, min_read_x) = params
    bam_path = folder + "/" + samfile
    sam_aln = pysam.AlignmentFile(bam_path, rb)
    barcodes_in_bam = {}
    for read_i in sam_aln.fetch(chr_name, start=start_x, stop=stop_x):
        if not read_i.has_tag(tag_name):
            continue

        barcode_i = read_i.get_tag(tag_name)
        if barcode_i in barcodes_in_bam:
            barcodes_in_bam[barcode_i] += 1
        else:
            barcodes_in_bam[barcode_i] = 1
    print("["+samfile+"] " + str(len(barcodes_in_bam)) + " barcodes were detected.")
    if barcodes == "NULL":
        print("["+samfile+"] Barcode list is not provided. Detetcting and filtering barcodes from the BAM file...")
        good_barcodes = [x for x in barcodes_in_bam if barcodes_in_bam[x] > min_read_x]
    elif barcodes != None:
        if not isinstance(barcodes, list):
            print("Error (Code 4): Parameter 'barcodes' should be 'NULL' or a list of 10X barcodes.")
            return -1
        print("["+samfile+"] Matching and filtering barcodes...")
        good_barcodes = [x for x in barcodes if x in barcodes_in_bam[x] and barcodes_in_bam[x] > min_read_x]
    else:
        print("Error (Code 3): Parameter 'barcodes' should be None of a list of 10X barcodes.")
        return -1
    print("["+samfile+"] " + str(len(good_barcodes)) + " barcodes were kept with reads count > " + str(min_read_x) + ".")
    return (samfile, good_barcodes, barcodes_in_bam)

def do_split_bam(params):
    (folder, samfile, rb, good_barcodes,tag_name,chr_name,threads) = params
    # print("[" + samfile + "]",len(good_barcodes))
    bam_path = folder + "/" + samfile
    print("["+samfile+"] Splitting BAM file...")

    suffix = "sam" if rb == "r" else "bam"

    output_folder = tempfile.mkdtemp()
    # output_folder = "example/split_bam/" + samfile[:-4]
    # os.mkdir(output_folder)

    chunk_size = int(500/threads)
    chunk_n = int(len(good_barcodes) / chunk_size)
    for chunk_i in range(chunk_n + 1):
        start_i = chunk_i * chunk_size
        end_i = (chunk_i+1) * chunk_size
        print("[" + samfile + "] "+str(chunk_i+1)+"/"+str(chunk_n+1))
        split_bam(bam_path, rb, good_barcodes[start_i:end_i], tag=tag_name, contigs=chr_name, output_tmp=output_folder)
        for barcode_i in good_barcodes[start_i:end_i]:
            pysam.index(output_folder + "/" + barcode_i + "." + suffix)

    print("[" + samfile + "] Splitting BAM file complete.")

    return (samfile,output_folder)


def read_sam(bam_folder, fmt="bam",sp="human",
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
    """
    :param folder: path(s) pointing to BAM alignment file(s).
    :param format: File type, bam or sam
    :param chr_name: Name of mitochondrial genome as specified in the BAM files.
    :param seq_type: Sequencing type:bulk or sc.
    :param combined_bam: If the BAM is merged file or not (It should be set to True for 10X or droplet data).
    :param tag_name: The name of the tag corresponding the cellular barcode. Default = "CB". For droplet scRNA-seq only.
    :param barcodes: The barcode list corresponding to the cells. For 10X genomics scRNA-seq data only.
    :param min_read: The minimum number of read counts to be considered a valid barcode (cell) in the analysis. Default = 1000. For droplet scRNAseq technologies only.
    :param max_depth: The maximum depth of reads considered at any position.
    :param min_baseq: The minimum read base quality below which the base is ignored.
    :param min_mapq: The minimum map quality below which the read is ignored.
    :return:
    """
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
    print(str(sam_n) + " BAM/SAM files were detected.")

    seq_len = mus_mt_len if sp.lower() =="mus" else human_mt_len

    ######
    if seq_type == "bulk":
        print("Each BAM file is considered as a sample and running on each BAM file separately...")
        print("Reading BAM/SAM files...")
        print(str(n_jobs) + ' processes are running ...')

        params = [(folder, sam_ls[i],i,sam_n, rb, chr_name, 0, seq_len, min_baseq, min_mapq) for i in range(sam_n)]
        pool = Pool(processes=n_jobs)
        results = pool.map(do_coverage_count, params)
        pool.close()

        samples = []
        for i in range(sam_n):
            sample_x = Sample(sam_ls[i][:-4], seq_type="bulk",sp=sp)
            sample_x.info = "Bulk sequencing data"
            (sample_x.A,sample_x.C,sample_x.G, sample_x.T, sample_x.coverage,sample_x.total_reads) = results[i]
            samples.append(sample_x)

        print("Cleaning memory...")
    elif seq_type == "sc":
        if not combined_bam:
            print("Each BAM file is considered as a cell sample and running on each BAM file separately...")

            now = datetime.now()  # current date and time
            date_time = now.strftime("%b/%d/%Y %a %H:%M:%S")
            print(date_time + "Reading BAM/SAM files...")
            print(str(n_jobs) + ' processes are running ...')

            params = [(folder, sam_ls[i],i,sam_n, rb, chr_name, 0, seq_len, min_baseq, min_mapq) for i in
                      range(sam_n)]
            pool = Pool(processes=n_jobs)
            results = pool.map(do_coverage_count, params)
            pool.close()

            samples = []
            cell_ls = []
            for i in range(sam_n):
                cell_x = Cell(sam_ls[i][:-4],sp=sp)
                (cell_x.A,cell_x.C,cell_x.G, cell_x.T, cell_x.coverage,cell_x.total_reads) = results[i]
                cell_ls.append(cell_x)

            sample_name = os.path.basename(folder)
            sample = Sample(sample_name, seq_type="sc",sp=sp)
            sample.info = """Single cell sequencing data.
            Each BAM file are considered as one cell sample.
            Sample name: """ + sample_name + """ (Folder name that contains the BAM files)
            Cell count: """ + str(sam_n)

            sample.cell_count = sam_n
            sample.cells = cell_ls

            samples.append(sample)
            print("Cleaning temporary files...")
        ## combined_bam
        else:
            print("The BAM file contains many cells, and barcodes will be used for cell count...")
            print(str(n_jobs) + ' processes are running...')

            now = datetime.now()  # current date and time
            date_time = now.strftime("%b/%d/%Y %a %H:%M:%S")
            print(date_time + " Reading BAM/SAM files...")

            samples = []

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
            params = [(folder, sam_ls[i], rb, good_barcodes[sam_ls[i]],tag_name,chr_name,n_jobs) for i in range(sam_n)]
            pool = Pool(processes=n_jobs)
            results = pool.map(do_split_bam, params)
            pool.close()
            for x_i in results:
                tmp_split_folder[x_i[0]] = x_i[1]

            ######
            for sam_i in sam_ls:
                cell_size = len(good_barcodes[sam_i])
                suffix = "sam" if rb == "r" else "bam"
                params = [(tmp_split_folder[sam_i], good_barcodes[sam_i][i] + "." + suffix, i,cell_size,
                           rb, chr_name, 0, seq_len, min_baseq, min_mapq) for i in range(cell_size)]
                pool = Pool(processes=n_jobs)
                results = pool.map(do_coverage_count, params)
                pool.close()

                cell_ls = []
                for i in range(cell_size):
                    cell_x = Cell(good_barcodes[sam_i][i],sp=sp)
                    (cell_x.A,cell_x.C,cell_x.G, cell_x.T, cell_x.coverage,cell_x.total_reads) = results[i]
                    cell_ls.append(cell_x)


                sample = Sample(sam_i[:-4], seq_type="sc",sp=sp)
                sample.info = """Single cell sequencing data. 
                Cells were stored in on BAM for each sample.
                Sample name: """ + sam_i[:-4] + """ (BAM/SAM file name)
                Barcodes count (After filtering/Total):""" +  str(cell_size) +"/"+ str(len(barcodes_in_bam[sam_i]))

                sample.cell_count = cell_size
                sample.cells = cell_ls

                samples.append(sample)

            print("Cleaning temporary files...")
            for name_i in tmp_split_folder:
                shutil.rmtree(tmp_split_folder[name_i])

    else:
        print("Error (Code 2): Parameter 'type' setting is wrong (Valid value: 'bulk' OR 'sc').")
        return -1

    print("Reading BAM/SAM files...Done!")
    if len(samples) == 1:
        return samples[0]
    else:
        return samples
