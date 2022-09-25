import matplotlib.pyplot as plt

import class_def
import plots
import processing


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
    for i, filename in enumerate(cv_result):
        CV = class_def.CyclicVoltammetry(filename, path_in, area)
        CV.solute_conc(vol_in, solute_conc)
        CV.peak_finder()
        CVList.append(CV)
        # plots.cv_plot(CV.cv['E'], CV.cv['j'], CV.filename, path_in)
        # plt.show()
    for cvs in CVList:  #implement further steps and plots
        print('Scan rate: {}, solute volume: {}'.format(cvs.rate, cvs.vol_sol))