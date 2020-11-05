import os
import pysam
try:
    import importlib.resources as pkg_resources
except ImportError:
    # Try backported to PY<37 `importlib_resources`.
    import importlib_resources as pkg_resources

from . import data

def read_fasta(file,sp="human"):
    tag = "HomoMT" if sp=="human" else "MusMT"
    fa_file = pysam.FastaFile(file)
    fa_seq=fa_file.fetch(tag)
    fa_file.close()
    return fa_seq

with pkg_resources.path(data, "human_mt.fasta") as f:
    data_path = f
human_mt_seq = read_fasta(data_path,sp="human")
human_mt_len = len(human_mt_seq)


with pkg_resources.path(data, "MusMT.fasta") as f:
    data_path = f
mus_mt_seq = read_fasta(data_path, sp="mus")
mus_mt_len = len(mus_mt_seq)
