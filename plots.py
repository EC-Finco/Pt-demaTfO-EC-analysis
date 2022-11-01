import os

import matplotlib
import matplotlib.pyplot as plt
import processing
import numpy as np


def cv_plot(pot, j, filename, path_in, export="n", peaks=None, peak_flag="no"):
    pathwd = path_in
    # plot cyclic voltammetry and peaks selected
    plt.title(filename)
    plt.plot(pot, j*1000, label="CV data")
    if peak_flag == 'yes':
        plt.plot(peaks.E_pc, peaks.j_pc*1000, "x", label="Cathodic peaks")
        plt.plot(peaks.E_pa, peaks.j_pa*1000, "x", label="Anodic peaks")
    plt.xlabel(r'$\mathit{E}$ / V vs. Pt', fontsize=20)
    plt.ylabel(r'$\mathit{j}$ / $mA \, cm^{-2}$', fontsize=20)
    plt.legend()
    plt.tight_layout()
    if export == 'y':
        filepath = path_in + '/preprocessed CVs/'
        os.chdir(filepath)
        plt.savefig(filename.replace('.txt', '.png'))
    os.chdir(pathwd)


def cv_volumes(cv_result, path_in, area, solute):  # plotting CVs with the same scan rate and different volume of solute
    pathexp = path_in + '/preprocessed CVs/'
    os.chdir(pathexp)
    plt.figure()
    plt.title('Voltammetries at 50 mV/s with different content of ' + solute)
    for i in cv_result:
        rate = processing.scanrate(i)
        if rate == 0.05:
            cv = processing.preproc(path_in, i, area)
            inversion, length_scans, scan_dir = processing.scan_wise(cv)
            peaks = processing.peak_finder(cv, inversion, length_scans, scan_dir, i)[0]
            half_wave_pot = peaks['E1/2'].mean()
            cv.E = cv.E - half_wave_pot
            plt.plot(cv.E, cv.j*1000, label=str(processing.volume_module(i)[0]))
    plt.xlabel(r'$\mathit{E}$ / V vs. $E_{1/2}$', fontsize=20)
    plt.ylabel(r'$\mathit{j}$ / mA cm$^{-2}$', fontsize=20)
    plt.legend(title=r'Solute volume ($\mathit{\mu L}$)')
    plt.tight_layout()
    plt.savefig('plots different solute volumes.png')
    plt.show()
    os.chdir(path_in)


def peak_cur_lin_sqrt(x, y_cat, y_an, fit_cat, fit_an, r2_cat, r2_an, figure_name, path_in):
    pathwd = path_in
    pathexp = pathwd + '/fittings/'
    os.chdir(pathexp)
    plt.figure()
    plt.plot(x, y_cat*1000, "x", label="cat peak")  # cathodic plot
    plt.plot(x, fit_cat.fittedvalues*1000, label="cat fit R$^{2}$ = %.2f" % r2_cat)
    plt.plot(x, y_an*1000, "x", label="an peak")  # anodic plot
    plt.plot(x, fit_an.fittedvalues*1000, label="an fit R$^{2}$ = %.2f" % r2_an)
    plt.xlabel("$v_{s}^{1/2}$ / [mV$^{1/2}$/s$^{1/2}$]", fontsize=20)
    plt.ylabel(r"$\mathit{j}$ / [mA cm$^{-2}$]", fontsize=20)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figure_name)
    plt.show()
    os.chdir(pathwd)


def second_plot(conc_cat, conc_an, regr_catRS, regr_anRS, slope_cat, slope_an, solute, path_in):
    # plot second regression
    pathexp = path_in + '/fittings/'
    os.chdir(pathexp)
    figurename = 'Second regression plot.png'
    plt.figure()
    plt.title('Final regression', fontsize=20)
    plt.plot(conc_cat*1000, slope_cat, "x", label="cat second fit, R$^{2}$ = %.2f" % regr_catRS.rsquared_adj)
    plt.plot(conc_cat*1000, regr_catRS.fittedvalues)
    plt.plot(conc_an*1000, slope_an, "x", label="an second fit, R$^{2}$ = %.2f" % regr_anRS.rsquared_adj)
    plt.plot(conc_an*1000, regr_anRS.fittedvalues)
    plt.xlabel('$C_{%s}$ / [M]' % solute, fontsize=20)
    plt.ylabel('$D^{1/2}C$ / [cm M s$^{-1/2}$]', fontsize=20)
    plt.legend()
    plt.tight_layout()
    plt.savefig(figurename)
    plt.show()
    os.chdir(path_in)


def lin_peak_cur_sqrt():
    plt.xlabel("$v_{s}^{1/2}$ / [mV$^{1/2}$/s$^{1/2}$]", fontsize=20)
    plt.ylabel(r"$\mathit{j_{p}}$ / [mA cm$^{-2}$]", fontsize=20)
    plt.legend()
    plt.tight_layout()


####### PLOTS FOR CHRONOAMPEROMETRIES ########
def ca_explore(CA):
    plt.figure()
    plt.plot(CA.data['Time'], label='Time')
    plt.plot(CA.data['Current Density'], label='Current')
    plt.plot(CA.data['Potential'], label='Potential')
    plt.xlabel('index')
    plt.tight_layout()
    plt.legend()
    plt.show()


def chrono_amp(CA, export):
    plt.figure()  # create figure with all CAs from one repetition
    plt.title(CA.filename)
    for u, data in enumerate(CA.CA):
        plt.plot(data['Time'], data['Current Density'], label='%.2f' % CA.U[u])
    plt.xlim(0.00, 0.01)
    plt.ylabel('$j$ / A cm$^{-2}$', fontsize=16)
    plt.xlabel('Time / s', fontsize=16)
    plt.legend()
    plt.tight_layout()
    if export == 'y':
        for ext in ['png', 'eps']:
            plotname = CA.filename[:-3] + ext
            plt.savefig(plotname)
    plt.show()


def c_diff(CA, export):
    plt.figure()
    plt.title('Integrated charge', fontsize=18)
    plt.plot(CA.U, CA.Cdiff)
    plt.ylabel('C$_{diff}$ / mF cm$^{-2}$', fontsize=16)
    plt.xlabel('E / V Vs. Pt', fontsize=16)
    plt.tight_layout()
    if export == 'y':
        for ext in ['png', 'eps']:
            plotname = 'C diff' + CA.filename[:-3] + ext
            plt.savefig(plotname)
            plt.show()


def chrono_coul(CA, export):
    plt.figure()  # create figure with all chrono coulometries from one repetition
    plt.title('Chronocoulometry')
    for u, data in enumerate(CA.CC):
        plt.plot(data, label='%.2f' % CA.U[u])
    plt.xlim(0.00, 0.02)
    plt.ylabel('$Q$ / C cm$^{-2}$')
    plt.xlabel('Time / s')
    plt.legend()
    plt.show()


def ca_3d(CA):
    # 3D plot
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    for u, data in enumerate(CA.CA):
        ax.plot_surface(data['Time'], data['Potential'], data['Current Density'], rstride=1, cstride=1,
                        cmap='viridis', edgecolor='none')


def fit_ca(t, y, y_fit, pot):
    plt.figure()
    plt.title('CA fitting at .2%f V' % pot)
    plt.plot(t, y, label="Data")
    plt.plot(t, y_fit, label="Fit")
    plt.legend()
    plt.show()
