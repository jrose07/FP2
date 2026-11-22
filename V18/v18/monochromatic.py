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
counts = np.genfromtxt("raw/Caesium.Spe", skip_header=12, max_rows=8192)
counts = counts / ((1*60*60) + (18*60) + 57)  # Messzeit Caesium: 1h18min57s
counts = counts - underground

channels = np.arange(len(counts))
energies = noms(m) * channels + noms(b)  # Channels in Energie
# Kanten der Kanäle für plt.stairs (ein Wert mehr als counts)
edges = noms(m) * (np.arange(len(counts) + 1) - 0.5) + noms(b)

plt.figure()
plt.stairs(counts, edges, fill=True, label="Cs-137")
plt.vlines(32,-10 , 0.1,color="red", linestyles="dashed",label="Ba")
plt.vlines(184,-10 , 0.1,color="red", linestyles="dashed",label="Rück")
plt.vlines(477,-10 , 0.1,color="red", linestyles="dashed",label="Com")
plt.vlines(662,-10 , 0.1,color="red", linestyles="dashed",label="Photopeak")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Count rate $A$ / s$^{-1}$")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "monochromatic_Cs.pdf")
plt.close()

# Für den Fit der Photolinie (Energie, Halbwertsbreite, Zehntelwertsbreite, Linieninhalt)
def gauss(K, h, mu, sigma, c):
    return h * np.exp(-(K - mu)**2 / (2 * sigma**2)) + c

# Photolinie: Gauß + Konstante in Kanälen fitten, danach in Energie umrechnen
K_max = np.argmax(counts[:8000])
window = 40  # Kanäle links und rechts vom Peak
K_fit = channels[K_max - window:K_max + window]
y_fit = counts[K_max - window:K_max + window]
p0 = [counts[K_max], K_max, 5, np.median(y_fit)]
p, cov_g = curve_fit(gauss, K_fit, y_fit, p0=p0)
h, mu, sigma, c = correlated_values(p, cov_g)
sigma = abs(sigma)

E_photo = m * mu + b
sigma_E = m * sigma
print(f"Photolinie: E = {E_photo} keV, sigma = {sigma_E} keV")

# Halbwerts- und Zehntelwertsbreite aus dem Gauß-Fit
FWHM_gauss = 2 * np.sqrt(2 * np.log(2)) * sigma_E
FWTM_gauss = 2 * np.sqrt(2 * np.log(10)) * sigma_E

# Halbwerts- und Zehntelwertsbreite direkt aus den Daten (linear interpoliert)
def width_at(fraction):
    level = noms(c) + fraction * noms(h)
    # links und rechts vom Maximum den ersten Kanal unter dem Level suchen
    left = K_max - np.argmax(counts[K_max::-1] < level)
    right = K_max + np.argmax(counts[K_max:] < level)
    K_left = np.interp(level, [counts[left], counts[left + 1]], [left, left + 1])
    K_right = np.interp(level, [counts[right], counts[right - 1]], [right, right - 1])
    return noms(m) * (K_right - K_left)

FWHM_data = width_at(1 / 2)
FWTM_data = width_at(1 / 10)
print(f"Halbwertsbreite:   Gauss = {FWHM_gauss} keV, Daten = {FWHM_data:.3f} keV")
print(f"Zehntelwertsbreite: Gauss = {FWTM_gauss} keV, Daten = {FWTM_data:.3f} keV")
print(f"Verhältnis FWTM/FWHM: Daten = {FWTM_data / FWHM_data:.3f}, Gauss = {np.sqrt(np.log(10) / np.log(2)):.3f}")

# Linieninhalt = Fläche unter dem Gauß (in Kanälen, daher Rate in 1/s)
N_photo = np.sqrt(2 * np.pi) * h * sigma
print(f"Linieninhalt Photolinie: N = {N_photo} 1/s")

K_fine = np.linspace(K_fit[0], K_fit[-1], 500)
plt.figure()
plt.plot(noms(m) * K_fit + noms(b), y_fit, ".", label="Data")
plt.plot(noms(m) * K_fine + noms(b), gauss(K_fine, *p), label="Gauss fit")
plt.axhline(noms(c) + noms(h) / 2, color="gray", ls="--", lw=0.8, label="Half maximum")
plt.axhline(noms(c) + noms(h) / 10, color="gray", ls=":", lw=0.8, label="Tenth maximum")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Count rate $A$ / s$^{-1}$")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "photopeak_Cs.pdf")
plt.close()
