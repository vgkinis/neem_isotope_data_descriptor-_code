import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, FixedFormatter
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import dates as mpl_dates
import string
import time
from datetime import datetime
import copy
import os
import subprocess
import sys
sys.path.append("/Users/vasilis")
import vaspy
from vaspy import syttensen
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# plt.rc('font', size = 12)
plt.rc('xtick', direction = 'in')
plt.rc('ytick', direction = 'in')
mylabelsize = 14
plt.rc('ytick', labelsize = mylabelsize)
plt.rc('xtick', labelsize = mylabelsize)
plt.rc('axes', labelsize = mylabelsize)
#comment

plt.close("all")
plt.ion()

all_stats = pd.read_excel("./neem_high_res_calcs/stats_all_neem_runs.xlsx", sheet_name = None, na_values = ["Nan"])
df_stats = all_stats["stats_NEEM_all_corr"].set_index(keys = 'time_stamp')
f = open("./neem_nd_all_corr_pangaea.txt", "r")
header_comments_pangaea_file = f.readlines()[0:92]
f.close()

f = open("./neem_nd_all_corr_pangaea_with_stats.txt", "w")
for jline in header_comments_pangaea_file:
    f.write(jline)
f.write("Depth\tage_GICC05\tage_GICC05_modelext\tage_AICC12\tGICC05_MCE\td18\tdD\t1-sigma_of_run_d18\t1-sigma_of_run_dD\taccuracy_check_offset_d18\taccuracy_check_offset_dD\n")
f.write("[m]\t[yb2k]\t[yb2k]\t[yb2k]\t[y]\t[permile]\t[permile]\t[permile]\t[permile]\t[permile]\t[permile]\n")


df_pangaea = pd.read_csv("./neem_nd_all_corr_pangaea.txt", sep = "\t", skiprows = 92).iloc[1:]
df_nd_all_cor = pd.read_csv("./neem_nd_all_corr.txt", header = 0, delimiter = "\t")
indexes_below_8k = np.arange(np.size(df_nd_all_cor["Depth"]))[32000:]

for i in indexes_below_8k:
    depth_i = df_nd_all_cor.iloc[i]["Depth"]
    depth_diff = np.around(depth_i - df_stats["top_depth"], 3)
    min_dif_positive = np.min(depth_diff[depth_diff>=0])
    j_on_stats = np.where(np.abs(depth_diff - min_dif_positive)<0.00001)[0][0]
    # print(depth_i, j_on_stats, df_stats.iloc[j_on_stats]["top_depth"])
    i_on_pangaea = np.where(np.abs(df_pangaea["Depth"].astype(float) - depth_i)<0.0001)
    if np.size(i_on_pangaea) == 0:
        continue
    data_row = "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.3f\t%0.2f\t%0.5f\t%0.5f\t%0.5f\t%0.5f" %(\
    df_pangaea["Depth"].to_numpy(dtype = float)[i_on_pangaea][0],\
    df_pangaea["age_GICC05"].to_numpy(dtype = float)[i_on_pangaea][0],\
    df_pangaea["age_GICC05_modelext"].to_numpy(dtype = float)[i_on_pangaea][0],\
    df_pangaea["age_AICC12"].to_numpy(dtype = float)[i_on_pangaea][0],\
    df_pangaea["GICC05_MCE"].to_numpy(dtype = float)[i_on_pangaea][0],\
    df_pangaea["d18"].to_numpy(dtype = float)[i_on_pangaea][0],\
    df_pangaea["dD"].to_numpy(dtype = float)[i_on_pangaea][0],\
    df_stats["std_std2_d18"].to_numpy(dtype = float)[j_on_stats],\
    df_stats["std_std2_dD"].to_numpy(dtype = float)[j_on_stats],\
    df_stats["offset_d18"].to_numpy(dtype = float)[j_on_stats],\
    df_stats["offset_dD"].to_numpy(dtype = float)[j_on_stats])
    if i%500==0:
        print(data_row)
    f.write(data_row)
    f.write("\n")

f.close()
