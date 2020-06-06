from scipy import signal
from scipy.fftpack import fft, fftshift
import numpy as np
from srxraylib.plot.gol import plot

def peak(x, c):
    return np.exp(-np.power(x - c, 2) / 16.0)

def lin_interp(x, y, i, half):
    return x[i] + (x[i+1] - x[i]) * ((half - y[i]) / (y[i+1] - y[i]))

def half_max_x(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[0], half),
            lin_interp(x, y, zero_crossings_i[1], half)]

def fwhm(x, y):
    half = max(y)/2.0
    signs = np.sign(np.add(y, -half))
    zero_crossings = (signs[0:-2] != signs[1:-1])
    zero_crossings_i = np.where(zero_crossings)[0]
    return [lin_interp(x, y, zero_crossings_i[1], half) - lin_interp(x, y, zero_crossings_i[0], half)]

en = (in_object_1.get_contents("xoppy_data")[:,0]).tolist()
ref = (in_object_1.get_contents("xoppy_data")[:,1]).tolist()
en_max = en[ref.index(max(ref))]
print("Max reflectivity: %f at: %d [eV]"%(ref[ref.index(max(ref))], en_max))

en_fwhm = fwhm(en, ref)[0]
print("FWHM: %f (%.2f %% BW)"%(en_fwhm, en_fwhm/en_max*100))
