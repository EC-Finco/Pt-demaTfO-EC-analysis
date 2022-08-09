import os

def folders_out(path_in, new_folders="n"):
    if new_folders == "y":
        path_smoothed = path_in + "/preprocessed CVs/"
        path_peaks = path_in + "/peak data/"
        os.mkdir(path_smoothed)
        os.mkdir(path_peaks)
    else:
        path_smoothed = path_in + "/preprocessed CVs/"
        path_peaks = path_in + "/peak data/"
    return path_smoothed, path_peaks

def exporter(cv_name, cv, peaks, path_in):
    path_smoothed, path_peaks = folders_out(path_in)
    file_cv = os.path.join(path_smoothed, cv_name)
    cv.to_csv(file_cv, sep="\t")
    file_peaks = os.path.join(path_peaks, cv_name)
    peaks.to_csv(file_peaks, sep="\t")
