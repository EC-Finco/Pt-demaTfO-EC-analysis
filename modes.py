import matplotlib.pyplot as plt
import pandas as pd

import plots
import processing, paths
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


def Randles_Sevcik(path_in, cv_result, export="n", vol_in="5", area="1"):
    solute_conc = float(input("Input solute concentration [mol/L]:\t"))
    conc, vol_solute, vol_solute_unit = processing.cv_features(cv_result, vol_in, solute_conc)
    conc_cat, conc_an, slope_cat, slope_an = processing.first_regression(cv_result, vol_solute_unit, path_in, export, area, conc)
    RS_df = processing.second_regression(conc_cat, conc_an, slope_cat, slope_an)
    print(RS_df)
    RS_df.to_csv('Randles-Sevcik analysis result.txt', sep='\t', header=True)


def cv(path_in, cv_result, export="n", area="1"):
    for i in cv_result:
        cv = processing.preproc(path_in, i, area)
        plots.cv_plot(cv.E, cv.j, i, export)
