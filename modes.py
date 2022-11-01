import matplotlib.pyplot as plt

import class_def
import plots
import processing
import numpy as np


def randles_sevcik(path_in, cv_result, export="n", vol_in="5", area="1"):
    solute_conc = float(input("Input solute concentration [mol/L]:\t"))
    solute = input('Type the solute species')
    conc, vol_solute, vol_solute_unit = processing.cv_features(cv_result, vol_in, solute_conc)
    # analysis
    conc_cat, conc_an, slope_cat, slope_an = processing.first_regression(cv_result, vol_solute_unit, path_in, export,
                                                                         area, conc)
    RS_df = processing.second_regression(conc_cat, conc_an, slope_cat, slope_an, path_in, solute)
    print(RS_df)
    RS_df.to_csv('Randles-Sevcik analysis result.txt', sep='\t', header=True)
    # plotting
    plots.cv_volumes(cv_result, path_in, area, solute)


def cvsurvey(path_in, cv_result, export="n", area="1"):
    for i in cv_result:
        cv = processing.preproc(path_in, i, area)
        plt.figure()
        plots.cv_plot(cv.E, cv.j, i, export)
        plt.show()


### ALTERNATIVE RANDLES-SEVCIK ANALYSIS
def randles_sevcik2(path_in, cv_result, export, vol_in, area):
    solute_conc = float(input("Input solute concentration [mol/L]:\t"))
    solute = input('Type the solute species')
    CVList = []
    CVVol = []
    LinListAn = []
    LinListCat = []
    for i, filename in enumerate(cv_result):
        CV = class_def.CyclicVoltammetry(filename, path_in, area)
        CV.solute_conc(vol_in, solute_conc)
        CV.peak_finder()
        CVList.append(CV)
        # plots.cv_plot(CV.cv['E'], CV.cv['j'], CV.filename, path_in)
        # plt.show()
    cond = processing.volumes_unique(CVList)
    # tuple of lists of unique values for volume (pos. 0) all values for respective quantities in pos 1
    for vol in cond[0]:  # collect sqrt(vs) and jp anodic and cathodic for each volume addition
        index = [idx for idx in range(len(CVList)) if cond[1][idx] == vol]
        for i in index:
            CVVol.append(CVList[i])
        rates, jpa, jpc = processing.jp_extraction(CVVol)
        LinAn = class_def.CurrentRateLin(rates, jpa, vol)
        LinCat = class_def.CurrentRateLin(rates, jpc, vol)
        LinListAn.append(LinAn)
        LinListCat.append(LinCat)
        CVVol = []  # fix order of data
    plt.figure()
    for LinCat in LinListCat:
        plt.plot(LinCat.rate_sq, LinCat.jp*1000, 'x', label=LinCat.vol)
        plt.plot(LinCat.rate_sq, LinCat.fitted*1000, label='%.2f' % LinCat.r2)
        plots.lin_peak_cur_sqrt()
    plt.show()
    # to complete the second regression and fix the plotting


# Chronoamperometry analysis
def chronoamp(path_in, ca_result, export, area):
    print('Found %d CAs in the folder' % len(ca_result))
    for i, filename in enumerate(ca_result):
        CA = class_def.ChronoAmperometry(filename, path_in, area)  # create object from a file
        plots.ca_explore(CA)
        export = 'y'
        plots.chrono_amp(CA, export)  # plots all the CAs performed at the different potentials
        plots.c_diff(CA, export)
        CA.chrono_fit()


