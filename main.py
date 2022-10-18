# -----------------------------------MODULE CHECKS------------------------------

# Check for modules, try to exit gracefully if not found
import sys
import importlib
import paths
import processing
import modes
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
import math
# Stop message from appearing
import warnings

warnings.filterwarnings("ignore", ".*GUI is implemented.*")
warnings.filterwarnings("ignore", ".*No labelled objects found.*")
# stable definition for code testing
new_folders = 'n'
mode = 'CA'
export = 'n'
diam = 8
# insert the path of the spectra
path_in = input("Type the path of files: ")
os.chdir(path_in)
# new_folders = input("create new folders for code output? [y/n]\t")
# mode = input("Select the kind of study that you want to perform: \n\t[CV]-simple CV plotting \n\t[P]-find peaks "
#            "\n\t[L]-linearize current/scan rate relation \n\t[RS]-solute concentration effect "
#            "\n\t[CA]-Chronoamperometry")
path_smoothed, path_peaks = paths.folders_out(path_in, new_folders)
# export = input("Export files? [y/n]\t")
# Find relevant CVs in folder
cv_result = [i for i in glob.glob('cv*.txt')]  # only case-insensitive in Windows
ca_result = [i for i in glob.glob('ca*.txt')]  # only case-insensitive in Windows
# diam = float(input("Input electrode diameter in millimeters:\t"))
Area = math.pi * (diam/20)**2  # area in cm^2
if mode == 'RS':
    vol_in = input("Insert the initial volume of the electrolyte in mL:\t")
    modes.randles_sevcik2(path_in, cv_result, export, vol_in, Area)
if mode == 'CV':
    modes.cvsurvey(path_in, cv_result, export, Area)
if mode == 'CA':
    modes.chronoamp(path_in, ca_result, export, Area)
