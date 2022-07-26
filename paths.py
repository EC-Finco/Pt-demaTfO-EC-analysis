import os


def folders_out(path_in, new_folders):
    if new_folders == "y":
        path_smoothed = path_in + "/preprocessed CVs/"
        path_peaks = path_in + "/peak data/"
        os.mkdir(path_smoothed)
        os.mkdir(path_peaks)
    else:
        path_smoothed = path_in + "/preprocessed spectra/"
        path_peaks = path_in + "/peak data/"
    return path_smoothed, path_peaks
