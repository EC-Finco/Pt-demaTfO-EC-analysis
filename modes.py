import matplotlib.pyplot as plt
import pandas as pd

import processing, paths
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score


def Randles_Sevcik(path_in, cv_result, export="n", vol_in="5", area="1"):
    solute_conc = float(input("Input solute concentration [mol/L]:\t"))
    conc, vol_solute, vol_solute_unit = processing.cv_features(cv_result, vol_in, solute_conc)
    conc_cat, conc_an, slope_cat, slope_an = processing.first_regression(cv_result, vol_solute_unit, path_in, export, area, conc)
    RS_array = processing.second_regression(conc_cat, conc_an, slope_cat, slope_an)
    RS_df = pd.DataFrame(RS_array)
    RS_df.columns = ['Slope', 'Intercept', 'R^2']
    RS_df.index = ['cathodic', 'anodic']
    print(RS_df)
