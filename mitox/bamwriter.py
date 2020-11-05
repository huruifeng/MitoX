import pysam

class BamWriter:
    def __init__(self, alignment, rb, barcodes, prefix):
        self.alignment = alignment
        self.rb=rb
        self.prefix = prefix
        self.barcodes = set(barcodes)
        self._out_files = {}

    def write_barcode_read(self, read, barcode):
        if barcode not in self.barcodes:
            return
        if barcode not in self._out_files:
            self._open_file_for_barcode(barcode)
        self._out_files[barcode].write(read)

    def _open_file_for_barcode(self, barcode):
        if self.rb=="r":
            suffix = "sam"
            w="w"
        else:
            suffix = "bam"
            w="wb"
        self._out_files[barcode] = pysam.AlignmentFile(
            f"{self.prefix}/{barcode}.{suffix}", "wb", template=self.alignment)
    def close_all(self):
        for barcode_i in self.barcodes:
            self._out_files[barcode_i].close()



def split_bam(input_bam,rb, barcodes, tag, contigs, output_tmp):
    """
     Split 10x barcoded BAM file into barcode-specific (cell specific) BAMs
    :param input_bam:  the 10x BAM file
    :param barcodes: a list of barcodes
    :param output_prefix: prefix of output files
    :return:
    """
    sam_aln = pysam.AlignmentFile(input_bam,rb)
    writer = BamWriter(alignment=sam_aln,rb=rb, barcodes=barcodes, prefix=output_tmp)

    for read in sam_aln.fetch(contigs):
        try:
            if not read.has_tag(tag):
                continue
            barcode = read.get_tag(tag)
            writer.write_barcode_read(read=read, barcode=barcode)
        except Exception as e:
            raise Exception("""Error (Code 6): Error in dealing with barcoded BAM file.""")

    writer.close_all()
