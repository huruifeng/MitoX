import os

import numpy as np
import pandas as pd

def mut_mean_std(x_df):
    res_df = pd.DataFrame(columns=["mean","std"], index=x_df.index)
    res_df['mean'] = x_df.mean(axis=1)
    res_df['std'] = x_df.std(axis=1)
    return res_df


def cal_dis_ij(af_ij,cov_ij,c=0):
    i = np.all((cov_ij > c), axis=0)
    is_all_zero = np.all((i == 0))
    if is_all_zero:
        d = 0.5
    else:
        af_ij = af_ij[:,i]
        d = np.mean(np.sqrt(np.absolute(af_ij[0]-af_ij[1])))
    return d


