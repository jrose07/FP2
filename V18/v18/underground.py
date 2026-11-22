import uncertainties.unumpy as unp
import numpy as np
from scipy.optimize import curve_fit
from addons import write, add, latex_float, tab_to_latex as tab2tex
from scipy.stats import linregress
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
import scipy.constants as const
import matplotlib.pyplot as plt
from scipy.signal import find_peaks


dir = "content/plots/"
dir_tab = "content/tables/"

# Spektrum laden: Zeile 13 bis 8204 der .Spe-Datei sind die 8192 Kanäle
counts = np.genfromtxt("raw/V18_Max_and_Julian/Untergrund.Spe", skip_header=12, max_rows=8192)
counts = counts / ((21*60*60) + (11*60))
channels = np.arange(len(counts))

def get_underground():
    return counts

peaks, props = find_peaks(counts, prominence=50)

print(peaks)
print(counts[peaks])

plt.figure()
plt.stairs(counts, np.arange(len(counts) + 1), fill=True, label="Eu-152")
plt.plot(peaks, counts[peaks], "x", color="red", label="peaks")
plt.xlabel("Channel")
plt.ylabel("Counts $N$")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "Underground.pdf")
plt.close()
