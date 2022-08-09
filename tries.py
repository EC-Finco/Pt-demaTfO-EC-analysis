import numpy as np

peak_cur_rate = np.zeros((5, 2, 2))
print(peak_cur_rate)
def peak_coll(cv, n_scans, scan_dir, inversion):
    scans_I = []  # data from each scan in a different array
    scans_E = []
    peaks_indexs = []
    peaks_I = []
    peaks_E = []
    peaks_pos = []
    for j in range(n_scans - 1):
        scans_I.append(np.array(0))
        scans_I[j] = cv.I[int(inversion[j]):int(inversion[j + 1])]
        scans_E.append(np.array(0))
        scans_E[j] = cv.E[int(inversion[j]):int(inversion[j + 1])]
        peaks_I.append(np.array(0))
        peaks_E.append(np.array(0))
        peaks_indexs.append(np.array(0))
    #  peak finder
    for i in range(n_scans - 1):
        if scan_dir[i] == 'anodic':
            peaks_pos = find_peaks(scans_I[i], distance=50)[0]
            scan = scans_I[i]
            peaks_I[i] = np.take(scan, peaks_pos)
            scan = scans_E[i]
            peaks_E[i] = np.take(scan, peaks_pos)
        elif scan_dir[i] == 'cathodic':
            peaks_pos = find_peaks(scans_I[i], distance=50)[0]
            scan = scans_I[i]
            peaks_I[i] = np.take(scan, peaks_pos)
            scan = scans_E[i]
            peaks_E[i] = np.take(scan, peaks_pos)
        peaks_indexs[i] = peaks_pos
        # to be continued...
        # peak selector
        # peaks_indexs[i] = peaks_indexs[i][scans_I[peaks_indexs[i]] == max(scans_I[peaks_indexs[i]])]
