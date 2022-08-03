import matplotlib.pyplot as plt

import processing, paths
import numpy as np


def Randles_Sevcik(path_in, cv_result, export="n", vol_in="5"):
    solute_conc = float(input("Input solute concentration"))
    conc, conc_list, vol_solute = processing.cv_features(cv_result, vol_in, solute_conc)
    peak_info = np.zeros((len(cv_result), 3))
    k = 0  # index of file being analyzed
    n_rates = []
    plot_test = "n"  # change to input to turn on cv plotting
    for j in vol_solute:
        j = str(j)
        j = j[:-2]
        l = 0
        for i in cv_result:
            if j in i:  # only runs in files with the same added volume
                file_input = i
                cv = processing.preproc(path_in, file_input)
                inversion, length_scans, scan_dir = processing.scan_wise(cv)
                peaks, peak_info[k, :] = processing.peak_info(cv, inversion, length_scans, scan_dir, file_input)
                k = k + 1
                l = l + 1
                if export == 'y':
                    paths.exporter(i, cv, peaks, path_in)
                if plot_test == "y":  # testing by plotting
                    plt.figure()
                    plt.plot(cv.E, cv.I)
                    plt.plot(peaks.E_pc, peaks.I_pc, "x", label="cat peak")
                    plt.plot(peaks.E_pa, peaks.I_pa, "x", label="an peak")
                    plt.show()
        n_rates.append(l)  # number of different scan rates probed for each addition
    plt.figure()
    plt.plot(peak_info[:, 0], peak_info[:, 1], "x", label="cat peak")
    plt.plot(peak_info[:, 0], peak_info[:, 2], "x", label="an peak")
    plt.show()
    # put plotting in outer for and add linear regression against v1/2
    # proceed with analysis