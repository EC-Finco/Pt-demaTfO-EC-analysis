# -----------------------------------MODULE CHECKS------------------------------

# Check for modules, try to exit gracefully if not found
import sys
import importlib
import paths
import processing
try:
    importlib.import_module('numpy')
    foundnp = True
except ImportError:
    foundnp = False
try:
    importlib.import_module('matplotlib')
    foundplot = True
except ImportError:
    foundplot = False
try:
    importlib.import_module('pandas')
    foundpd = True
except ImportError:
    foundplot = False
if not foundnp:
    print("Numpy is required. Exiting")
    sys.exit()
if not foundplot:
    print("Matplotlib is required. Exiting")
    sys.exit()
if not foundpd:
    print("Pandas is required. Exiting")
    sys.exit()
try:
    importlib.import_module('pybaselines')
    foundbase = True
except ImportError:
    foundbase = False
if not foundbase:
    print("Pybaselines is required. Exiting")
    sys.exit()
# -------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pybaselines
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks, peak_prominences, peak_widths
# -------------------------------------------------------------------------------

import os
import glob

# Stop message from appearing
import warnings

warnings.filterwarnings("ignore", ".*GUI is implemented.*")
warnings.filterwarnings("ignore", ".*No labelled objects found.*")


# insert the path of the spectra
path_in = input("Type the path of CV: ")
os.chdir(path_in)
new_folders = input("create new folders for code output? [y/n]")
path_smoothed, path_peaks = paths.folders_out(path_in, new_folders)

# Find relevant CVs in folder
cv_result = [i for i in glob.glob('cv*.txt')]  # only case-insensitive in Windows
file_input = cv_result[0]
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
# plt.plot(cv.E, cv.I)
# plt.plot(cv.E, cv.dI_dE)
# plt.plot(cv.dE)
# plt.show()
# find inversions
inversion = np.zeros(len(cv.dE))  # indexing variable
for i in range(len(cv.dE)-1):
    if cv.dE[i]*cv.dE[i+1] <= 0:
        inversion[i] = i
inversion = inversion[(inversion > 0)]  # removes empty elements
# split scans
n_scans = len(inversion)  # counts scans
length_scans = np.zeros(n_scans)  # measures length of scan
scan_dir = []  # label scan direction
for i in range(n_scans-1):
    length_scans[i] = int(inversion[i+1]-inversion[i])
    if cv.dE[inversion[i]+1] >= 0:
        scan_dir.append('anodic')
    else:
        scan_dir.append('cathodic')
an_peaks = pd.DataFrame()
cat_peaks = pd.DataFrame()
for j in range(n_scans - 1):
    index_start = int(inversion[j])
    index_end = int(inversion[j+1])
    scan = cv.iloc[index_start:index_end]
    index = []
    for i in range(int(length_scans[j])):
        index.append(str(i+1))
    scan.index = index
    if scan_dir[j] == 'anodic':
        peak_index = np.argmax(scan.I)
        print(peak_index)
        print(scan)
        an_peaks.at[j,'Index'] = peak_index
        an_peaks.at[j,'E_p'] = scan.E[peak_index]
        an_peaks.at[j,'I_p'] = scan.I[peak_index]
    elif scan_dir[j] = 'cathodic':
        peak_index = np.argmax(scan.E)  # second derivative
print(an_peaks)
# # print(length_scans, scan_dir)  # use to probe the output
# scan = np.array(1)
# an_peaks = np.array(1)  # column 0 is index, column 1 is peak pot and column 2 is peak current
# cat_peaks = np.array(1)
# for j in range(n_scans - 1):
#     scan[:, 0] = cv.E[inversion[j]:inversion[j+1]]
#     scan[:, 1] = cv.I[inversion[j]:inversion[j+1]]
#     if scan_dir[j] == 'anodic':
#         an_peaks[j, 0] = np.argmax(scan[:, 1])
#         an_peaks[j, 1:2] = scan[an_peaks[j, 0], :]
#     # elif scan_dir[j] == 'cathodic':


        # baselining to remove discharge or indexing maxima



