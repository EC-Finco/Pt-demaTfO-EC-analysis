import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import statsmodels.api as smapi


def cv_plot(pot, j, filename, export="n", peaks=None, peak_flag="no"):
    # plot cyclic voltammetry and peaks selected
    plt.figure()
    plt.title(filename)
    plt.plot(pot, j*1000, label="CV data")
    if peak_flag == 'yes':
        plt.plot(peaks.E_pc, peaks.j_pc*1000, "x", label="Cathodic peaks")
        plt.plot(peaks.E_pa, peaks.j_pa*1000, "x", label="Anodic peaks")
    plt.xlabel('$\mathit{E}$ / V vs. Pt', fontsize=20)
    plt.ylabel('$\mathit{j}$ / $mA \, cm^{-2}$', fontsize=20)
    plt.legend()
    plt.tight_layout()
    if export == 'y':
        plt.savefig(filename.replace('.txt', '.png'))
    plt.show()


def peak_cur_lin_sqrt(x, y_cat, y_an, fit_cat, fit_an, r2_cat, r2_an, figure_name):
    plt.figure()
    plt.plot(x, y_cat*1000, "x", label="cat peak")  # cathodic plot
    plt.plot(x, fit_cat.fittedvalues*1000, label="cat fit R$^{2}$ = %.2f" % r2_cat)
    plt.plot(x, y_an*1000, "x", label="an peak")  # anodic plot
    plt.plot(x, fit_an.fittedvalues*1000, label="an fit R$^{2}$ = %.2f" % r2_an)
    plt.xlabel("$v_{s}^{1/2}$ / [mV$^{1/2}$/s$^{1/2}$]", fontsize=20)
    plt.ylabel("$\mathit{j}$ / [mA cm$^{-2}$]", fontsize=20)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()


def second_plot(conc_cat, conc_an, regr_catRS, regr_anRS, slope_cat, slope_an):
    # plot second regression
    figurename = 'Second regression plot.png'
    plt.figure()
    plt.title('Final regression', fontsize=20)
    plt.plot(conc_cat*1000, slope_cat, "x", label="cat second fit, R$^{2}$ = %.2f" % regr_catRS.rsquared_adj)
    plt.plot(conc_cat*1000, regr_catRS.fittedvalues)
    plt.plot(conc_an*1000, slope_an, "x", label="an second fit, R$^{2}$ = %.2f" % regr_anRS.rsquared_adj)
    plt.plot(conc_an*1000, regr_anRS.fittedvalues)
    plt.xlabel('$C_{HClO_{4}}$ / [M]', fontsize=20)
    plt.ylabel('$D^{1/2}C$ / [cm M s$^{-1/2}$]', fontsize=20)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figurename)
    plt.show()
