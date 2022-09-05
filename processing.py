import os.path
import numpy as np
import pandas as pd
import glob
import re
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy.ndimage import uniform_filter1d
import math
import paths
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import statsmodels.api as smapi

import plots

matplotlib.rc('text', usetex = True)

def preproc(path_in, file_input, area):
    path_file = path_in + "/" + file_input
    cv = pd.read_csv(path_file, delimiter="	", usecols=['Potential applied (V)', 'WE(1).Current (A)']
                     , engine='python', skipfooter=1)
    cv.columns = ["E", "I"]
    cv['j'] = uniform_filter1d(cv.I, 10) / area  # A cm-2
    # calculating derivatives
    dE = np.diff(cv.E)
    dj = np.diff(cv.j)
    dj = uniform_filter1d(dj, 10)
    dj_dE = dj / dE
    # moving to dataframe
    dE = np.insert(dE, 0, 0)
    dj = np.insert(dj, 0, 0)
    dj_dE = np.insert(dj_dE, 0, 0)
    cv['dE'] = dE
    cv['dj'] = dj
    cv['dj_dE'] = dj_dE
    return cv


def scan_wise(cv):
    # find inversions
    inversion = np.zeros(len(cv.dE))  # indexing variable
    for i in range(len(cv.dE) - 1):
        if cv.dE[i] * cv.dE[i + 1] <= 0:
            inversion[i] = i
    inversion = inversion[(inversion > 0)]  # removes empty elements
    # split scans
    n_scans = len(inversion)  # counts scans
    length_scans = np.zeros(n_scans)  # measures length of scan
    scan_dir = []  # label scan direction
    for i in range(n_scans - 1):
        length_scans[i] = int(inversion[i + 1] - inversion[i])
        if cv.dE[inversion[i] + 1] >= 0:
            scan_dir.append('anodic')
        else:
            scan_dir.append('cathodic')
    return inversion, length_scans, scan_dir


def peak_finder(cv, inversion, length_scans, scan_dir, file_input):
    n_scans = len(inversion)
    peaks = pd.DataFrame()
    cycle = 0
    for j in range(n_scans - 2):
        j = j+1
        index_start = int(inversion[j])
        index_end = int(inversion[j + 1])
        scan = cv.iloc[index_start:index_end]
        index = []
        for i in range(int(length_scans[j])):  # generating the index list to reindex each scan
            index.append(str(i + 1))
        scan.index = index
        if scan_dir[j] == 'anodic':  # anodic peak finder
            peak_index = np.argmax(scan.j)
            peaks.at[cycle, 'index_a'] = peak_index
            peaks.at[cycle, 'E_pa'] = scan.E[peak_index]
            peaks.at[cycle, 'j_pa'] = scan.j[peak_index]
        elif scan_dir[j] == 'cathodic':  # cathodic peak finder
            peaklist, _ = find_peaks(-scan.j)
            peak_indexes = np.array(peaklist, dtype=int)
            peak_potentials = scan.E[peak_indexes]
            if j >= 1:
                peak_indexes = peak_indexes[peak_potentials < peaks['E_pa'][cycle]]
            if peak_indexes.size != 0:
                peak_index = peak_indexes[0]
                peaks.at[cycle, 'index_c'] = peak_index
                peaks.at[cycle, 'E_pc'] = scan.E[peak_index]
                peaks.at[cycle, 'j_pc'] = scan.j[peak_index]
            cycle = cycle + 1
    peaks.at[:, 'DE_p'] = peaks['E_pa'][:] - peaks['E_pc'][:]
    peaks.at[:, 'Vol_sol'] = volume_module(file_input)[0]
    peaks.at[:, 'Scan_rate'] = scanrate(file_input)
    peak_info = np.zeros(4)
    peak_info[0] = peaks.Scan_rate[0]  # working attribution with first cycle
    peak_info[1] = math.sqrt(peaks.Scan_rate[0])
    peak_info[2] = peaks.j_pa[0]
    peak_info[3] = peaks.j_pc[0]
    return peaks, peak_info


# RANDLES SEVCIK SPECIFIC


#  procedure to list the features as solute concentration and scan rate from filenames
def volume_module(filename):
    regex = re.compile(r' \d+ uL')
    regex2 = re.compile(r' \d+uL')
    vol_unit = str()
    x = regex.findall(filename)  # string with uL
    if not x:
        x = regex2.findall(filename)  # alternative string with uL but no space between number and unit
        y = str(x[0])
        vol_unit = y
        y = y[:-2]
    else:
        y = str(x[0])
        vol_unit = y
        y = y[:-3]
    vol_sol = int(y)
    return vol_sol, vol_unit


def cv_features(cv_result, vol_ini, solute_conc):  # fix this one
    vol_in = float(vol_ini)
    vol_array = np.zeros(len(cv_result))
    vol_unit = []
    vol_unit_list = []
    j = 0
    y = str()
    unit = str()
    for i in cv_result:
        vol_array[j] = volume_module(i)[0]
        vol_unit.append(volume_module(i)[1])
        j = j + 1
    vol_solute_list = []
    for i in vol_array:  # removes duplicates
        if i not in vol_solute_list:
            vol_solute_list.append(i)
    for i in vol_unit:  # removes duplicates
        if i not in vol_unit_list:
            vol_unit_list.append(i)
    vol_solute = np.array(vol_solute_list)
    tot_vol = (vol_solute / 1000) + vol_in  # solute volume is in uL, solution in mL
    n_mol = (vol_solute / 1000000) * solute_conc
    conc = n_mol / tot_vol  # conc in mol/mL
    return conc, vol_solute_list, vol_unit_list


def scanrate(filename):
    regex = re.compile(r'\d+ mV[ ]{0,1}s-1')
    x2 = regex.findall(filename)  # string with uL
    y2 = str(x2[0])
    y2 = y2[:-7]
    rate = int(y2) / 1000  # conversion to V/s
    return rate


def first_regression(cv_result, vol_solute_list, path_in, export, Area, conc):
    peak_info = np.zeros((len(cv_result), 4))  # array containing one sampled peak for each file
    k = 0  # index of file being analyzed
    n_vol = 0  # index of the volume being analyzed
    slope_cat = np.zeros((len(vol_solute_list)))
    slope_an = np.zeros((len(vol_solute_list)))
    r2_cat = np.zeros((len(vol_solute_list)))
    r2_an = np.zeros((len(vol_solute_list)))
    plot_test = "y"  # change to turn on cv plotting
    for j in vol_solute_list:
        j = str(j)
        l = 0  # index of files with same solute concentration to be analyzed
        for i in cv_result:
            if j in i:  # only runs in files with the same added volume
                file_input = i
                cv = preproc(path_in, file_input, Area)
                inversion, length_scans, scan_dir = scan_wise(cv)
                peaks, peak_info[k, :] = peak_finder(cv, inversion, length_scans, scan_dir, file_input)
                k = k + 1
                l = l + 1
                if export == 'y':
                    paths.exporter(i, cv, peaks, path_in)
                if plot_test == "y":  # testing by plotting
                    plots.cv_plot(cv.E, cv.j, file_input, export, peaks, peak_flag="yes")
        peak_same_vol = peak_info[k - l:k, :]  # sampling from the general file
# first linear regression
        x = peak_same_vol[:, 1]
        y_cat = peak_same_vol[:, 2]
        y_an = peak_same_vol[:, 3]
        # reshaping
        # x = x1.reshape(-1, 1)
        # Create linear regression object
        fit_cat = gen_regr(x, y_cat)
        fit_an = gen_regr(x, y_an)
        # calculate R squared for both
        r2_cat[n_vol] = fit_cat.rsquared_adj
        r2_an[n_vol] = fit_an.rsquared_adj
        # extract fitting parameters
        slope_cat[n_vol] = fit_cat.params[1]
        slope_an[n_vol] = fit_an.params[1]
        # plot
        figure_name = 'Current linearization - ' + str(j) + ' uL.png'
        plots.peak_cur_lin_sqrt(x, y_cat, y_an, fit_cat, fit_an, r2_cat[n_vol], r2_an[n_vol], figure_name)
        n_vol = n_vol + 1
        k = 0
    # remove data with low R^2
    conc_cat1 = conc[r2_cat > 0.9]
    conc_an1 = conc[r2_an > 0.9]
    slope_cat = slope_cat[r2_cat > 0.9]
    slope_an = slope_an[r2_an > 0.9]
    conc_cat = conc_cat1.reshape(-1, 1)
    conc_an = conc_an1.reshape(-1, 1)
    return conc_cat, conc_an, slope_cat, slope_an


def second_regression(conc_cat, conc_an, slope_cat, slope_an):
    slope_cat = slope_cat * ((8.314 * 298) ** (1 / 2)) / (0.4463 * (96485 ** (1.5)))
    slope_an = slope_an * ((8.314 * 298) ** (1 / 2)) / (0.4463 * (96485 ** (1.5)))
    regr_catRS = gen_regr(conc_cat, slope_cat)
    regr_anRS = gen_regr(conc_an, slope_an)
    plots.second_plot(conc_cat, conc_an, regr_catRS, regr_anRS, slope_cat, slope_an)
    Dval_cat = (regr_catRS.params[1] ** 2)
    Derr_cat = abs(2 * regr_catRS.params[1] * regr_catRS.bse[1])
    C_0val_cat = abs(regr_catRS.params[0] / regr_catRS.params[1])
    C_0err_cat = abs(C_0val_cat * ((regr_catRS.bse[1] / regr_catRS.params[1]) + (regr_catRS.bse[0] / regr_catRS.params[0])))
    Dval_an = abs(regr_anRS.params[1] ** 2)
    Derr_an = abs(2 * regr_anRS.params[1] * regr_anRS.bse[1])
    C_0val_an = abs(- regr_anRS.params[0] / regr_anRS.params[1])
    C_0err_an = abs(C_0val_an * ((regr_anRS.bse[1] / regr_anRS.params[1]) + (regr_anRS.bse[0] / regr_anRS.params[0])))
    columns = ['D Value', 'D Error', 'C_0 Value', 'C_0 Error', 'Adj. R sq.']
    index = ['Cathodic', 'Anodic']
    values = [[Dval_cat, Derr_cat, C_0val_cat, C_0err_cat, regr_catRS.rsquared_adj],
              [Dval_an, Derr_an, C_0val_an, C_0err_an, regr_anRS.rsquared_adj]]
    RS_df = pd.DataFrame(values, columns=columns, index=index)
    return RS_df


def gen_regr(x, y):  # linear regression with statsmodel, add fitted points
    x = smapi.add_constant(x)  # adding the intercept term
    res_fit = smapi.OLS(y, x).fit()
    return res_fit



