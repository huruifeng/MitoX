
from .global_data import human_mt_seq
from .global_data import human_mt_len

from .global_data import mus_mt_seq
from .global_data import mus_mt_len

import pandas as pd
import numpy as np

import anndata as ad

np.random.seed(12)


class Sample:
    def __init__(self, name,sp="human",seq_type=None):
        self.species = sp
        self.name = name
        # A,T,C,G.N rate along the sequence
        self.A = []
        self.T = []
        self.C = []
        self.G = []
        self.N = []
        self.coverage = []
        self.total_reads=0

        self.mtDNA = ""
        self.mt_len = 0

        if self.species.upper() == "HUMAN":
            self.mtDNA = human_mt_seq
            self.mt_len = human_mt_len
        else:
            self.mtDNA=mus_mt_seq
            self.mt_len=mus_mt_len

        self.seq_type = seq_type
        self.cell_count = 0
        self.cells = []
        self.info = "NULL"

    def get_name(self):
        return self.name

    def get_ref_base(self,pos):
        pos = pos -1
        return self.mtDNA[pos]

    def get_coverage(self, pos):
        pos = pos - 1
        if self.seq_type == "bulk":
            return self.coverage[pos]
        elif self.seq_type == "sc":
            print("Sample contains several cells, please use sample[s_i].cells[c_i].get_coverage(pos=pos_i) to get a position-specific coverage.")

    def get_allele_count(self,pos):
        pos = pos - 1
        if self.seq_type == "bulk":
            ref_i = self.mtDNA[pos]
            return {"ref": ref_i.upper(), "A": self.A[pos], "C": self.C[pos], "G": self.G[pos], "T": self.T[pos]}
        elif self.seq_type == "sc":
            print("Sample contains several cells, please use sample[s_i].cells[c_i].get_allele_count(pos=pos_i) to get a position-specific allele counts.")

    def get_allele_rate(self, pos):
        pos = pos - 1
        if self.seq_type == "bulk":
            ref_i = self.mtDNA[pos]
            covegage_i = self.coverage[pos]
            return {"ref": ref_i.upper(),
                    "A": float(self.A[pos])/covegage_i,
                    "C": float(self.C[pos])/covegage_i,
                    "G": float(self.G[pos])/covegage_i,
                    "T": float(self.T[pos])/covegage_i}
        elif self.seq_type == "sc":
            print(
                "Sample contains several cells, please use sample[s_i].cells[c_i].get_allele_rate(pos=pos_i) to get a position-specific allele rates.")

    def gen_mut_profile(self):
        mut_dict = {}
        mut_cov = {}
        index_ls = []
        index_pos = []
        base = ["A","C","G","T"]
        for i in range(self.mt_len):
            for base_i in base:
                if self.mtDNA[i]!=base_i:
                    index_ls.append(str(i+1)+self.mtDNA[i]+">"+base_i)
                    index_pos.append(str(i+1)+self.mtDNA[i]+">"+base_i)
        if self.seq_type == "sc":
            for cell_i in self.cells:
                mut_dict[cell_i.name] = []
                mut_cov[cell_i.name] = []
                for i in range(self.mt_len):
                    ref_i = self.mtDNA[i]
                    if ref_i == "A":
                        AC_rate = float(cell_i.C[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        AG_rate = float(cell_i.G[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        AT_rate = float(cell_i.T[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        mut_dict[cell_i.name] += [AC_rate,AG_rate,AT_rate]
                    elif ref_i == "C":
                        CA_rate = float(cell_i.A[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        CG_rate = float(cell_i.G[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        CT_rate = float(cell_i.T[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        mut_dict[cell_i.name] += [CA_rate, CG_rate, CT_rate]
                    elif ref_i == "G":
                        GA_rate = float(cell_i.A[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        GC_rate = float(cell_i.C[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        GT_rate = float(cell_i.T[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        mut_dict[cell_i.name] += [GA_rate, GC_rate, GT_rate]
                    elif ref_i == "T":
                        TA_rate = float(cell_i.A[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        TC_rate = float(cell_i.C[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        TG_rate = float(cell_i.G[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        mut_dict[cell_i.name] += [TA_rate, TC_rate, TG_rate]
                    else:
                        NA_rate = float(cell_i.A[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        NC_rate = float(cell_i.C[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        NG_rate = float(cell_i.G[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        NT_rate = float(cell_i.T[i]) / cell_i.coverage[i] if cell_i.coverage[i] != 0 else 0.0
                        mut_dict[cell_i.name] += [NA_rate,NC_rate, NG_rate, NT_rate]

                    if ref_i == "N":
                        mut_cov[cell_i.name] += [cell_i.coverage[i], cell_i.coverage[i], cell_i.coverage[i],cell_i.coverage[i]]
                    else:
                        mut_cov[cell_i.name] += [cell_i.coverage[i], cell_i.coverage[i], cell_i.coverage[i]]
        else:
            print("""the data is bulk sequencing data, please use the function in gen_mut_profile(samples) in mitox.util. e.g.
                     from mitox.utils import gen_mut_profile
                     mut_profile = gen_mut_profile(sample list)""")

        profile_df = pd.DataFrame(data=mut_dict, index=index_ls)
        df_coverage = pd.DataFrame(data=mut_cov, index=index_pos)

        adata = ad.AnnData(profile_df.T)
        adata.varm["coverage"] = df_coverage
        adata.uns["species"] = self.species.lower()
        adata.uns['mut'] = True

        return adata

class Cell:
    def __init__(self, name,sp="human"):
        self.species = sp
        self.name = name
        # A,T,C,G.N rate along the sequence
        self.A = []
        self.T = []
        self.C = []
        self.G = []
        self.N = []
        self.coverage = []
        self.total_reads = 0

        self.mtDNA = ""
        self.mt_len = 0

        if self.species.upper() == "HUMAN":
            self.mtDNA = human_mt_seq
            self.mt_len = human_mt_len
        else:
            self.mtDNA = mus_mt_seq
            self.mt_len = mus_mt_len

    def get_name(self):
        return self.name

    def get_ref_base(self,pos):
        pos = pos - 1
        return self.mtDNA[pos]

    def get_coverage(self, pos):
        pos = pos - 1
        return self.coverage[pos]

    def get_allele_count(self, pos):
        pos = pos - 1
        ref_i = self.mtDNA[pos]
        return {"ref": ref_i.upper(), "A": self.A[pos], "C": self.C[pos], "G": self.G[pos], "T": self.T[pos]}

    def get_allele_rate(self, pos):
        pos = pos - 1
        ref_i = self.mtDNA[pos]
        covegage_i = self.coverage[pos]
        return {"ref": ref_i.upper(),
                "A": float(self.A[pos]) / covegage_i,
                "C": float(self.C[pos]) / covegage_i,
                "G": float(self.G[pos]) / covegage_i,
                "T": float(self.T[pos]) / covegage_i}


