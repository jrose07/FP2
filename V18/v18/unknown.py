import uncertainties.unumpy as unp
import numpy as np
from scipy.optimize import curve_fit
from addons import write, add, latex_float, tab_to_latex as tab2tex
from scipy.stats import linregress
from uncertainties import ufloat, correlated_values
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
import scipy.constants as const
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from underground import get_underground

underground = get_underground()
dir = "content/plots/"
dir_tab = "content/tables/"

# Energiekalibrierung aus calibration.py: E = m * K + b
m = ufloat(0.11950, 0.00005)  # keV/Kanal
b = ufloat(-0.98, 0.23)       # keV

# Spektrum laden: Zeile 13 bis 8204 der .Spe-Datei sind die 8192 Kanäle
counts = np.genfromtxt("raw/V18_Max_and_Julian/was_das_question.Spe", skip_header=12, max_rows=8192)
counts = counts / ((1*60*60) + (2*60) + 15)  # Messzeit Barium: 1h02min15s
counts = counts - underground

channels = np.arange(len(counts))
energies = noms(m) * channels + noms(b)  # Channels in Energie
# Kanten der Kanäle für plt.stairs (ein Wert mehr als counts)
edges = noms(m) * (np.arange(len(counts) + 1) - 0.5) + noms(b)

peaks, props = find_peaks(counts[:8000], prominence=0.05, width=4)
print("Peaks:",peaks)

plt.figure()
plt.stairs(counts, edges, fill=True, label="Ba-133")
plt.plot(energies[peaks], counts[peaks], "x", color="red", label="peaks")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Count rate $A$ / s$^{-1}$")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "unknown.pdf")
plt.close()
