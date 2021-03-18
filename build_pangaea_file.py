import numpy as np
import pandas as pd
import openpyxl
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, FixedFormatter
from matplotlib.backends.backend_pdf import PdfPages
import string
import time
import copy
import os
import subprocess
plt.rc('text', usetex=True)
plt.rc('font', family='serif')
# plt.rc('font', sans-serif = ['Helvetica'])
plt.rc('xtick', direction = 'in')
plt.rc('ytick', direction = 'in')
plt.rc('font', size = 12)
plt.rc('lines', linewidth = 1)
plt.close("all")
plt.ion()

class PangaeaBuilder():
    def __init__(self):
        return

    def read_iso_data(self):
        self.data_iso = np.genfromtxt("./neem_nd_all_corr.txt", skip_header = 1, usecols = (3,4,5))
        print(self.data_iso)
        print(np.shape(self.data_iso))
        self.z_iso_all = self.data_iso[:, 0]
        self.d18_iso_all = self.data_iso[:, 1]
        self.dD_iso_all = self.data_iso[:, 2]
        index_b8 = np.where(self.z_iso_all>1210.45)[0][0]
        self.z_iso_b8 = self.z_iso_all[index_b8:]
        self.d18_iso_b8 = self.d18_iso_all[index_b8:]
        self.dD_iso_b8 = self.dD_iso_all[index_b8:]
        return

    def plot_depth_all(self):
        f0, axess = plt.subplots(nrows = 1, ncols = 1 , num = 200, figsize = (8,6), tight_layout = True)
        ax0 = axess
        ax0.plot(self.z_iso_all, self.d18_iso_all, linewidth = 0.5, color = "b")
        ax0.plot(self.z_iso_b8, self.d18_iso_b8, linewidth = 0.5, color = "r")
        ax0.set_ylabel(r"$\delta^{18}\mathrm{O}$ [\textperthousand]")
        ax0.set_xlabel("Depth [m]")
        return

    def read_age_data(self):
        f = open("./NEEM_GICC05modelext_AICC12_5cm.txt", "r")
        self.age_file_comments_header = f.readlines()[0:82]
        f.close()
        print(self.age_file_comments_header)
        self.data_age = np.genfromtxt("./NEEM_GICC05modelext_AICC12_5cm.txt", skip_header = 84)
        print(self.data_age)
        self.z_age_all = self.data_age[:, 0]
        self.age_GICC05 = self.data_age[:, 1]
        self.age_GICC05_modelext = self.data_age[:, 2]
        self.age_AICC12 = self.data_age[:, 3]
        self.age_GICC05_MCE = self.data_age[:, 4]
        index_b8 = np.where(self.z_age_all>1210.45)[0][0]
        self.z_age_b8 = self.z_age_all[index_b8:]
        self.age_GICC05_b8 = self.age_GICC05[index_b8:]
        self.age_GICC05_modelext_b8 = self.age_GICC05_modelext[index_b8:]
        self.age_AICC12_b8 = self.age_AICC12[index_b8:]
        self.age_GICC05_MCE_b8 = self.age_GICC05_MCE[index_b8:]
        return

    def build_pangaea_file(self):
        print(np.shape(self.z_iso_b8))
        print(np.shape(self.z_age_b8))
        nans_to_append = np.zeros(self.z_iso_b8.size - self.z_age_b8.size) -999.99
        print(nans_to_append)
        self.age_GICC05_b8_append = np.append(self.age_GICC05_b8, nans_to_append)
        self.age_GICC05_modelext_b8_append = np.append(self.age_GICC05_modelext_b8, nans_to_append)
        self.age_AICC12_b8_append = np.append(self.age_AICC12_b8, nans_to_append)
        self.age_GICC05_MCE_b8_append = np.append(self.age_GICC05_MCE_b8 , nans_to_append)

        dataout = np.transpose(np.vstack((self.z_iso_b8, self.age_GICC05_b8_append, \
        self.age_GICC05_modelext_b8_append, self.age_AICC12_b8_append, \
        self.age_GICC05_MCE_b8_append, self.d18_iso_b8, self.dD_iso_b8)))
        print(dataout)
        f = open("./neem_nd_all_corr_pangaea.txt", "w")
        f.write("NEEM ice core stable water isotope record in the interval 1210.5 m (8000 years b2k) - 2536.5 m.\n")
        f.write("The record is measured on discrete samples with a depth resolution of 0.05 m\n")
        f.write("using Cavity Ring Down Spectroscopy. Measurements are reported in\n\
permile on the SMOW-SLAP scale using a two-fixed-point calibration.\n\
Columns 1-sigma_of_run_d18 and 1-sigma_of_run_dD contain the precision of \n\
each individual analysis run the data-point belongs to.\n\
Columns accuracy_check_offset_d18 and accuracy_check_offset_dD contain the \n\
accuracy check metrics of each analysis run.\n\
**All depth, age and isotope values refer to the bottom of each interval**")
        f.write("\n" + 40*"-" + "\n")
        for j in self.age_file_comments_header:
            f.write(j)
        f.write("Depth\tage_GICC05\tage_GICC05_modelext\tage_AICC12\tGICC05_MCE\td18\tdD\n")
        f.write("[m]\t[yb2k]\t[yb2k]\t[yb2k]\t[y]\t[permile]\t[permile]\n")
        np.savetxt(f, dataout, fmt = "%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.2f\t%0.3f\t%0.2f")
        f.close()
        return

    def plot_age_data(self):
        f0, axess = plt.subplots(nrows = 1, ncols = 1 , num = 300, figsize = (8,6), tight_layout = True)
        ax0 = axess
        ax0.scatter(self.age_GICC05_b8_append[np.where(self.age_GICC05_b8_append!=-999.99)], \
        self.d18_iso_b8[np.where(self.age_GICC05_b8_append!=-999.99)], s = 0.5, color = "b", label = "GICC05")
        ax0.scatter(self.age_GICC05_modelext_b8_append[np.where(self.age_GICC05_modelext_b8_append!=-999.99)], \
        self.d18_iso_b8[np.where(self.age_GICC05_modelext_b8_append!=-999.99)], color = "r", s = 0.5, label = "GICC05-modelext")
        ax0.scatter(self.age_AICC12_b8_append[np.where(self.age_AICC12_b8_append!=-999.99)], \
        self.d18_iso_b8[np.where(self.age_AICC12_b8_append!=-999.99)], color = "g", s = 0.5, label = "AICC12")
        ax0.set_ylabel(r"$\delta^{18}\mathrm{O}$ [\textperthousand]")
        ax0.set_xlabel("Age [yb2k]")
        ax0.legend()
        return



pb = PangaeaBuilder()
pb.read_iso_data()
pb.plot_depth_all()
pb.read_age_data()
pb.build_pangaea_file()
pb.plot_age_data()
