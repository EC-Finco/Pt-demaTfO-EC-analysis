import numpy as np
import pandas as pd
import glob
import re
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy.ndimage import uniform_filter1d


def preproc(path_in, file_input):
    path_file = path_in + "/" + file_input
    cv = pd.read_csv(path_file, delimiter="	", header=1, names=["E", "I"], engine='python', skipfooter=1)
    cv.I = uniform_filter1d(cv.I, 10)
    # calculating derivatives
    dE = np.diff(cv.E)
    dI = np.diff(cv.I)
    dI = uniform_filter1d(dI, 10)
    dI_dE = dI / dE
    # moving to dataframe
    dE = np.insert(dE, 0, 0)
    dI = np.insert(dI, 0, 0)
    dI_dE = np.insert(dI_dE, 0, 0)
    cv['dE'] = dE
    cv['dI'] = dI
    cv['dI_dE'] = dI_dE
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


def peak_info(cv, inversion, length_scans, scan_dir, file_input):
    n_scans = len(inversion)
    peaks = pd.DataFrame()
    cycle = 0
    for j in range(n_scans - 1):
        index_start = int(inversion[j])
        index_end = int(inversion[j + 1])
        scan = cv.iloc[index_start:index_end]
        index = []
        for i in range(int(length_scans[j])):  # generating the index list to reindex each scan
            index.append(str(i + 1))
        scan.index = index
        if scan_dir[j] == 'anodic':  # anodic peak finder
            peak_index = np.argmax(scan.I)
            peaks.at[cycle, 'Index_a'] = peak_index
            peaks.at[cycle, 'E_pa'] = scan.E[peak_index]
            peaks.at[cycle, 'I_pa'] = scan.I[peak_index]
        elif scan_dir[j] == 'cathodic':  # cathodic peak finder
            peaklist, _ = find_peaks(-scan.I)
            peak_indexes = np.array(peaklist, dtype=int)
            peak_potentials = scan.E[peak_indexes]
            if j >= 1:
                peak_indexes = peak_indexes[peak_potentials < peaks['E_pa'][cycle]]
            if peak_indexes.size != 0:
                peak_index = peak_indexes[0]
                peaks.at[cycle, 'Index_c'] = peak_index
                peaks.at[cycle, 'E_pc'] = scan.E[peak_index]
                peaks.at[cycle, 'I_pc'] = scan.I[peak_index]
            cycle = cycle + 1
    peaks.at[:, 'DE_p'] = peaks['E_pa'][:] - peaks['E_pc'][:]
    peaks.at[:, 'Vol_sol'] = volume_append(file_input)
    peaks.at[:, 'Scan_rate'] = scanrate(file_input)
    peak_info = np.zeros((3))
    peak_info[0] = peaks.Scan_rate[0]  # working attribution
    peak_info[1] = peaks.I_pa[0]
    peak_info[2] = peaks.I_pc[0]
    return peaks, peak_info

# RANDLES SEVCIK SPECIFIC

#  procedure to list the features as solute concentration and scan rate from filenames
def cv_features( cv_result, vol_ini, solute_conc):
    regex = re.compile(r'\d+ uL')
    vol_in = float(vol_ini)
    vol_array = np.zeros(len(cv_result))
    j = 0
    y = str()
    for i in cv_result:
        x = regex.findall(i)  # string with uL
        y = str(x[0])
        y = y[:-3]
        vol_array[j] = int(y)
        j = j+1
    vol_solute_list = []
    for i in vol_array:  #removes duplicates
        if i not in vol_solute_list:
            vol_solute_list.append(i)
    vol_solute = np.array(vol_solute_list)
    tot_vol = (vol_solute/1000) + vol_in  # solute volume is in uL
    n_mol = (vol_solute/1000000) * solute_conc
    conc = n_mol / (tot_vol/1000)
    conc_list = (vol_array/1000000) * solute_conc / ((vol_array/1000) + vol_in/1000)
    return conc, conc_list, vol_solute_list


def scanrate(filename):
    regex = re.compile(r'\d+ mV s-1')
    rate = 0
    y2 = str()
    x2 = regex.findall(filename)  # string with uL
    y2 = str(x2[0])
    y2 = y2[:-7]
    rate = int(y2)
    return rate




def volume_append(file):
    regex = re.compile(r'\d+ uL')
    y = str()
    x = regex.findall(file)  # string with uL
    y = str(x[0])
    y = y[:-3]
    vol_sol = int(y)
    return vol_sol