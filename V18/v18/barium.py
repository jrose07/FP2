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
from calibration import params_Q, cov_Q, m, b

underground = get_underground()
dir = "content/plots/"
dir_tab = "content/tables/"

# Energiekalibrierung E = m * K + b kommt aus calibration.py

# Spektrum laden: Zeile 13 bis 8204 der .Spe-Datei sind die 8192 Kanäle
counts = np.genfromtxt("raw/V18_Max_and_Julian/Barium.Spe", skip_header=12, max_rows=8192)
counts = counts / ((1*60*60) + (4*60) + 6)  # Messzeit Barium: 1h04min06s
counts = counts - underground

channels = np.arange(len(counts))
energies = noms(m) * channels + noms(b)  # Channels in Energie
# Kanten der Kanäle für plt.stairs (ein Wert mehr als counts)
edges = noms(m) * (np.arange(len(counts) + 1) - 0.5) + noms(b)

peaks, props = find_peaks(counts[:8000], prominence=0.03, width=4)



plt.figure()
plt.stairs(counts, edges, fill=True, label="Ba-133")
plt.plot(energies[peaks], counts[peaks], "x", color="red", label="peaks")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Count rate $A$ / s$^{-1}$")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "barium.pdf")
plt.close()

# Nur Linien über 150 keV (Q-Fit gilt erst ab dort)
peaks_fit = peaks[energies[peaks] > 150]

def gauss(K, h, mu, sigma, c):
    return h * np.exp(-(K - mu)**2 / (2 * sigma**2)) + c

window = 40  # Kanäle links und rechts vom Peak
N_line = []
E_peak = []
fig, axs = plt.subplots(1, len(peaks_fit), figsize=(3 * len(peaks_fit), 3.5))
for K, ax in zip(peaks_fit, axs):
    K_fit = channels[K - window:K + window]
    y_fit = counts[K - window:K + window]
    p0 = [counts[K], K, 5, np.median(y_fit)]
    p, cov_g = curve_fit(gauss, K_fit, y_fit, p0=p0)
    h, mu, sigma, c = correlated_values(p, cov_g)

    # Fläche unter dem Gauß = Linieninhalt (Rate in 1/s, weil counts schon durch t geteilt ist)
    N_line.append(np.sqrt(2 * np.pi) * h * abs(sigma))
    E_peak.append(m * mu + b)

    K_fine = np.linspace(K_fit[0], K_fit[-1], 300)
    ax.plot(energies[K_fit], y_fit, ".", markersize=3, label="Data")
    ax.plot(noms(m) * K_fine + noms(b), gauss(K_fine, *p), label="Gauss fit")
    ax.set_title(f"{noms(m * mu + b):.1f} keV")
    ax.set_xlabel("Energy $E$ / keV")
axs[0].set_ylabel("Count rate $A$ / s$^{-1}$")
axs[0].legend()
fig.tight_layout()
fig.savefig(dir + "barium_gauss_fits.pdf")
plt.close(fig)

N_line = np.array(N_line)
E_peak = np.array(E_peak)
for E, N in zip(E_peak, N_line):
    print(f"E = {E} keV   N = {N} 1/s")

# Emissionswahrscheinlichkeiten der Ba-133 Linien (NNDC), gleiche Reihenfolge wie peaks_fit
E_lit = np.array([276.3989, 302.8508, 356.0129, 383.8485])  # keV
W = np.array([7.16, 18.34, 62.05, 8.94]) / 100

# Q(E) = q_a * E^q_b direkt aus calibration.py (Parameter mit Kovarianzmatrix)
q_a, q_b = correlated_values(params_Q, cov_Q)
Q = np.array([q_a * E**q_b for E in E_peak])

# Raumwinkel mit dem Abstand der Barium-Quelle (gleich gemessen wie bei Europium)
a = 7.9 - 4.665 + 1.5  # Abstand Quelle - Detektor in cm
r = 4.5 / 2            # Radius des Ge-Kristalls in cm
Om = 0.5 * (1 - (a / (np.sqrt(a**2 + r**2))))  # Omega / 4pi
print("Omega/4pi =", Om)

# Aktivität A = N / (Q * W * Omega/4pi)
A = N_line / (Q * W * Om)
for E, q, act in zip(E_lit, Q, A):
    print(f"E = {E:7.2f} keV   Q = {q}   A = {act} Bq")

A_mean = np.mean(noms(A))
A_mean_err = np.std(noms(A), ddof=1) / np.sqrt(len(A))
print(f"Mittlere Aktivität: A = {A_mean:.0f} +/- {A_mean_err:.0f} Bq")
