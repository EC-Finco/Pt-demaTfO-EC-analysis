# -----------------------------------MODULE CHECKS------------------------------

# Check for modules, try to exit gracefully if not found
import sys
import importlib
import paths
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
path_in = input("Type the path of spectra: ")
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
plt.plot(cv.dE)
plt.show()
#split scans
inversion = [0 for i in range(len(cv.dE) - 1)]
for i in range(len(cv.dE)-1):
    sgn = cv.dE[i] * cv.dE[i+1]
    if sgn < 0:
        inversion[i] = i
inversion = inversion[inversion != 0]
print(dE)