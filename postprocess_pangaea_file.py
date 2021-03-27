import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, FixedLocator, FixedFormatter
from matplotlib.backends.backend_pdf import PdfPages
plt.rc('text', usetex=True)
plt.rc('font', family='sans-serif', size = 16)
plt.rc('text.latex', preamble=r'\usepackage{cmbright}')
plt.rc('xtick', direction = 'in')
plt.rc('ytick', direction = 'in')
mylabelsize = 16
plt.rc('ytick', labelsize = mylabelsize)
plt.rc('xtick', labelsize = mylabelsize)
plt.rc('axes', labelsize = mylabelsize)

plt.close("all")
plt.ion()


def simplify_pangaea(pangaea_fname, plotit = False):
    """
    Converts the PANGAEA file to one of a simpler form for loading/reading
    output file is ./simplify_pangaea.out
    """
    datain = np.genfromtxt(pangaea_fname, skip_header = 45, delimiter = "\t", filling_values = -999.99)
    datain
    print(datain)
    comment_header = """#NEEM ice core high resolution water isotope record
#Simplified version of the dataset https://doi.org/10.1594/PANGAEA.925552
#published in DOI:"""
    new_header = "\nDepth\tage_GICC05\tage_GICC05_modelext\tage_AICC12\tGICC05_MCE\td18\tdD\t1-sigma_of_run_d18\t1-sigma_of_run_dD\taccuracy_check_offset_d18\taccuracy_check_offset_dD\n"
    f = open("./simplify_pangaea.out", "w")
    f.write(comment_header)
    f.write(new_header)
    fmt_out = "%0.2f\t%0.6f\t%0.6f\t%0.6f\t%0.6f\t%0.3f\t%0.2f\t%0.5f\t%0.5f\t%0.5f\t%0.5f"
    np.savetxt(f, datain, fmt = fmt_out)
    f.close()

    if plotit:
        f0, axess = plt.subplots(nrows = 1, ncols = 1 , num = 510, figsize = (10,6), tight_layout = True)
        f0.suptitle("NEEM 0.05 m resolution record", size = 18)
        ax0 = axess
        ax0.plot(datain[:,0], datain[:,6], linewidth = 1, color = "r")
        ax0.set_ylabel(r'$\delta\mathrm{D}$ [\textperthousand]')
        ax0.set_xlabel(r'Depth [m]')

    return

def interpolate_on_age(pangaea_fname, plotit = False, interp_interval_ky = 5e-3):
    """
    Generates files with d18 and dD signals interpolated on an equal step
    time scale interpolated according to variable interp_interval_ky
    Generates ./interpolate_on_age_GICC05.out and ./interpolate_on_age_GICC05_AICC12.out
    """
    datain = np.genfromtxt(pangaea_fname, skip_header = 45, delimiter = "\t", filling_values = -999.99)
    depth = datain[:,0]
    age_GICC05 = datain[:,1]
    age_GICC05_modelext = datain[:,2]
    age_AICC12 = datain[:,3]
    mce_GICC05 = datain[:,4]
    d18 = datain[:,5]
    dD = datain[:,6]

    indexes_GICC05 = np.where(age_GICC05!=-999.99)[0]
    indexes_GICC05_modelext = np.where(age_GICC05_modelext!=-999.99)[0]
    indexes_AICC12 = np.where(age_AICC12!=-999.99)[0]
    indexes_GICC05_total = np.concatenate((indexes_GICC05, indexes_GICC05_modelext))
    indexes_AICC12_total = np.concatenate((indexes_GICC05, indexes_AICC12))

    age_GICC05_total = np.concatenate((age_GICC05[indexes_GICC05], age_GICC05_modelext[indexes_GICC05_modelext]))
    d18_GICC05_total = d18[indexes_GICC05_total]
    dD_GICC05_total = dD[indexes_GICC05_total]
    age_AICC12_total = np.concatenate((age_GICC05[indexes_GICC05], age_AICC12[indexes_AICC12]))
    d18_AICC12_total = d18[indexes_AICC12_total]
    dD_AICC12_total = dD[indexes_AICC12_total]


    indexes_sort_GICC05_total = np.argsort(age_GICC05_total)
    indexes_sort_AICC12_total = np.argsort(age_AICC12_total)

    age_GICC05_interp = np.arange(np.ceil(age_GICC05_total[indexes_sort_GICC05_total][0]*1000)/1000, \
        age_GICC05_total[indexes_sort_GICC05_total][-1], interp_interval_ky)
    d18_GICC05_interp = np.interp(age_GICC05_interp, age_GICC05_total[indexes_sort_GICC05_total], \
        d18_GICC05_total[indexes_sort_GICC05_total])
    dD_GICC05_interp = np.interp(age_GICC05_interp, age_GICC05_total[indexes_sort_GICC05_total], \
        dD_GICC05_total[indexes_sort_GICC05_total])

    age_AICC12_interp = np.arange(np.ceil(age_AICC12_total[indexes_sort_AICC12_total][0]*1000)/1000, \
        age_AICC12_total[indexes_sort_AICC12_total][-1], interp_interval_ky)
    d18_AICC12_interp = np.interp(age_AICC12_interp, age_AICC12_total[indexes_sort_AICC12_total], \
        d18_AICC12_total[indexes_sort_AICC12_total])
    dD_AICC12_interp = np.interp(age_AICC12_interp, age_AICC12_total[indexes_sort_AICC12_total], \
        dD_AICC12_total[indexes_sort_AICC12_total])

    if plotit:
        f0, axess = plt.subplots(nrows = 1, ncols = 1 , num = 510, figsize = (12,6), tight_layout = True)
        f0.suptitle("NEEM high resolution record vs age", size = 18)
        ax0 = axess
        ax0.plot(age_GICC05_total[indexes_sort_GICC05_total], dD_GICC05_total[indexes_sort_GICC05_total], \
            linewidth = 0.5, color = "indianred", label = "GICC05+GICC05ext")
        ax0.plot(age_AICC12_total[indexes_sort_AICC12_total], dD_AICC12_total[indexes_sort_AICC12_total], \
            linewidth = 0.5, color = "royalblue", label = "GICC05+AICC12")
        ax0.plot(age_GICC05_interp, dD_GICC05_interp, linewidth = 1, color = "mediumorchid", \
            label = "GICC05+GICC06ext interp: %0.2e ky" %interp_interval_ky)
        ax0.plot(age_AICC12_interp, dD_AICC12_interp, linewidth = 1, color = "forestgreen", \
            label = "GICC05+AICC12 interp: %0.2e ky" %interp_interval_ky)
        ax0.set_ylabel(r'$\delta\mathrm{D}$ [\textperthousand]')
        ax0.set_xlabel(r'Age [ky b2k]')
        ax0.legend(frameon = False, fontsize = 10)

    comment_header = """#NEEM ice core high resolution water isotope record
#the dataset https://doi.org/10.1594/PANGAEA.925552 interpolated on time
#published in DOI: """
    fout_GICC05 = open("./interpolate_on_age_GICC05.out", "w")
    dataout = np.transpose(np.vstack((age_GICC05_interp, d18_GICC05_interp, dD_GICC05_interp)))
    fout_GICC05.write("Age_GICC05+GICC05ext\td18\tdD\n")
    np.savetxt(fout_GICC05, dataout, fmt = "%0.6f\t%0.3f\t%0.2f")
    fout_GICC05.close()

    fout_AICC12 = open("./interpolate_on_age_GICC05_AICC12.out", "w")
    dataout = np.transpose(np.vstack((age_AICC12_interp, d18_AICC12_interp, dD_AICC12_interp)))
    fout_AICC12.write("Age_GICC05+AICC12\td18\tdD\n")
    np.savetxt(fout_AICC12, dataout, fmt = "%0.6f\t%0.3f\t%0.2f")
    fout_AICC12.close()



    return

# simplify_pangaea("./Gkinis-etal_2020.tab", plotit = False)
interpolate_on_age("./Gkinis-etal_2020.tab", plotit = False, interp_interval_ky = 2e-2)
