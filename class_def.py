import os.path
import numpy as np
import pandas as pd
import glob
import re
from scipy.signal import find_peaks, peak_prominences, peak_widths
from scipy import integrate
from scipy.ndimage import uniform_filter1d
import math
import paths
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import statsmodels.api as smapi
from scipy.optimize import curve_fit
import plots


class CyclicVoltammetry:
    def __init__(self, filename, path_in, area):
        self.peaks = None
        self.conc = None
        self.vol_sol = None
        self.vol_unit = None
        self.filename = filename
        self.path_file = path_in + "/" + filename
        self.cv = pd.read_csv(self.path_file, delimiter="	", usecols=['Potential applied (V)', 'WE(1).Current (A)'],
                              engine='python', skipfooter=1)
        self.cv.columns = ["E", "I"]
        self.cv['j'] = uniform_filter1d(self.cv.I, 10) / area  # A cm-2
        # calculating derivatives
        dE = np.diff(self.cv.E)
        dj = np.diff(self.cv.j)
        dj = uniform_filter1d(dj, 10)
        dj_dE = dj / dE
        # moving to dataframe
        dE = np.insert(dE, 0, 0)
        dj = np.insert(dj, 0, 0)
        dj_dE = np.insert(dj_dE, 0, 0)
        self.cv['dE'] = dE
        self.cv['dj'] = dj
        self.cv['dj_dE'] = dj_dE
        # find inversions
        self.inversion = np.zeros(len(self.cv.dE))  # indexing variable
        for i in range(len(self.cv.dE) - 1):
            if self.cv.dE[i] * self.cv.dE[i + 1] <= 0:
                self.inversion[i] = i
        self.inversion = self.inversion[(self.inversion > 0)]  # removes empty elements
        # split scans
        self.n_scans = len(self.inversion)  # counts scans
        self.length_scans = np.zeros(self.n_scans)  # measures length of scan
        self.scan_dir = []  # label scan direction
        for i in range(self.n_scans - 1):
            self.length_scans[i] = int(self.inversion[i + 1] - self.inversion[i])
            if self.cv.dE[self.inversion[i] + 1] >= 0:
                self.scan_dir.append('anodic')
            else:
                self.scan_dir.append('cathodic')
        regex = re.compile(r'\d+ mV[ ]{0,1}s-1')
        x2 = regex.findall(filename)  # string with uL
        y2 = str(x2[0])
        y2 = y2[:-6]
        self.rate = int(y2) / 1000  # conversion to V/s

    def solute_conc(self, vol_in, solute_conc):
        regex = re.compile(r' \d+ uL')
        regex2 = re.compile(r' \d+uL')
        x = regex.findall(self.filename)  # string with uL
        if not x:
            x = regex2.findall(self.filename)  # alternative string with uL but no space between number and unit
            y = str(x[0])
            self.vol_unit = y[-2:]
            y = y[:-2]
        else:
            y = str(x[0])
            self.vol_unit = y[-3:]
            y = y[:-3]
        self.vol_sol = int(y)
        tot_vol = (self.vol_sol / 1000) + float(vol_in)  # solute volume is in uL, solution in mL
        n_mol = (self.vol_sol / 1000000) * float(solute_conc)
        self.conc = n_mol / tot_vol  # conc in mol/mL

    def peak_finder(self):
        self.n_scans = len(self.inversion)
        self.peaks = pd.DataFrame()
        cycle = 0  # counts in which cycle the CV is being analyzed
        for j in range(self.n_scans - 1):
            j = j
            index_start = int(self.inversion[j])
            index_end = int(self.inversion[j + 1])
            scan = self.cv.iloc[index_start:index_end]
            index = []
            for i in range(int(self.length_scans[j])):  # generating the index list to reindex each scan
                index.append(str(i + 1))
            scan.index = index
            if self.scan_dir[j] == 'anodic':  # anodic peak finder
                peak_index = np.argmax(scan.j)
                self.peaks.at[cycle, 'index_a'] = peak_index
                self.peaks.at[cycle, 'E_pa'] = scan.E[peak_index]
                self.peaks.at[cycle, 'j_pa'] = scan.j[peak_index]
            elif self.scan_dir[j] == 'cathodic':  # cathodic peak finder
                peak_list, _ = find_peaks(-scan.j)
                peak_indexes = np.array(peak_list, dtype=int)
                peak_potentials = scan.E[peak_indexes]
                if j >= 1:
                    peak_indexes = peak_indexes[peak_potentials < self.peaks['E_pa'][cycle]]
                if peak_indexes.size != 0:
                    peak_index = peak_indexes[0]
                    self.peaks.at[cycle, 'index_c'] = peak_index
                    self.peaks.at[cycle, 'E_pc'] = scan.E[peak_index]
                    self.peaks.at[cycle, 'j_pc'] = scan.j[peak_index]
                cycle = cycle + 1
        self.peaks.at[:, 'DE_p'] = self.peaks['E_pa'][:] - self.peaks['E_pc'][:]
        self.peaks.at[:, 'E1/2'] = (self.peaks['E_pa'][:] + self.peaks['E_pc'][:]) / 2


class ChronoAmperometry:
    def __init__(self, filename, path_in, area):
        self.Cdiff = None
        self.Charge = []
        self.CA = []
        self.CC = []
        self.U = []
        self.filename = filename
        self.path_file = path_in + "/" + filename
        self.data_raw = pd.read_csv(self.path_file,
                                    delimiter="	", na_values=['ND'], engine='python', header=0).dropna()
        self.data = self.data_raw[['Corrected time (s)', 'WE(1).Current (A)', 'WE(1).Potential (V)']].copy()
        self.data.columns = ['Time', 'Current', 'Potential']
        self.data['Time'].astype(float)
        self.data['Current'].astype(float)
        self.data['Potential'].astype(float)
        self.data['Current'] = self.data['Current'] / area
        self.data.columns = ['Time', 'Current Density', 'Potential']
        self.restart_idx = self.data.index[self.data['Time'] == 0].tolist()
        # self.restart_idx = self.restart_idx[:-1]  # remove last CA which is usually faulty
        self.restart_idx = [x for i, x in enumerate(self.restart_idx) if i != 19]
        # removed 20th element since on that day it was wrong
        for j, i in enumerate(self.restart_idx):
            self.CA.append(self.data.loc[i:i+1199, 'Time':'Potential'])  # slicing the raw data into the list of
            # dataframes
            self.U.append(self.data.loc[i+1199, 'Potential'])
            # rewrite the index
            self.CA[j].index = list(range(1200))
        self.U = np.array(self.U)
        self.tare()
        self.integration()
        
    def tare(self):
        min_j = np.zeros(len(self.restart_idx))
        for i, ca in enumerate(self.CA):
            min_j[i] = ca['Current Density'].iloc[-1]
            self.CA[i]['Current Density'] = self.CA[i]['Current Density'] - min_j[i]

    # compute integral
    def integration(self):
        for ca in self.CA:
            self.CC.append(integrate.cumtrapz(ca['Current Density'], ca['Time'], initial=0))
        for i, cc in enumerate(self.CC):
            self.Charge.append(cc[-1])
        self.Charge = np.array(self.Charge, dtype=float)
        self.Cdiff = self.Charge * 1000 / 0.025  # potential step 25 mV, calculates differential capacitance in mF

    def chrono_fit(self):  # fitting only on the first 0,01 seconds
        for ca in self.CA:
            print(ca[1:119])
            pot = np.average(ca.loc[:119, 'Potential'])
            print(pot)
            t = np.array(ca.loc[:119, 'Time'])
            y = np.array(ca.loc[:119, 'Current Density'])
            param, param_cov = curve_fit(decay_func, t, y)
            y_fit = decay_func(t, param[0], param[1], param[2], param[3])
            plots.fit_ca(t, y, y_fit, pot)


class CurrentRateLin:  # input arrays containing data from CVs with same solute concentration
    def __init__(self, rate, jp, vol):
        self.vol = vol
        self.jp = jp
        self.rate_sq = np.sqrt(rate)
        x = self.rate_sq
        y = self.jp
        # Create linear regression object
        x = smapi.add_constant(x)  # adding the intercept term
        self.fit = smapi.OLS(y, x).fit()
        # calculate R squared for both
        self.r2 = self.fit.rsquared_adj
        # extract fitting parameters
        self.slope = self.fit.params[1]
        self.fitted = self.fit.fittedvalues


class DiffStudy:
    def __init__(self, slope, conc):
        self.slope = slope
        self.conc = conc
        y = self.slope * ((8.314 * 298) ** (1 / 2)) / (0.4463 * (96485 ** 1.5))
        x = self.conc
        # Create linear regression object
        x = smapi.add_constant(x)  # adding the intercept term
        fit = smapi.OLS(y, x).fit()
        # calculate R squared for both
        self.D_val = (fit.params[1] ** 2)
        self.D_err = abs(2 * fit.params[1] * fit.bse[1])
        self.C_0val = abs(fit.params[0] / fit.params[1])
        self.C_0err = abs(self.C_0val * ((fit.bse[1] / fit.params[1]) + (fit.bse[0] / fit.params[0])))


def decay_func(t, R_EDL, C_EDL, P1, P2):
    PhiS = 0.025  # set amplitude of step
    return (PhiS / R_EDL) * np.exp(- t / (R_EDL * C_EDL)) + P1 * np.exp(- P2 * t)
