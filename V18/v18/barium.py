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
plt.ylabel("Counts per second / s$^{-1}$")
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
for K in peaks_fit:
    K_fit = channels[K - window:K + window]
    y_fit = counts[K - window:K + window]
    p0 = [counts[K], K, 5, np.median(y_fit)]
    p, cov_g = curve_fit(gauss, K_fit, y_fit, p0=p0)
    h, mu, sigma, c = correlated_values(p, cov_g)

    # Fläche unter dem Gauß = Linieninhalt (Rate in 1/s, weil counts schon durch t geteilt ist)
    N_line.append(np.sqrt(2 * np.pi) * h * abs(sigma))
    E_peak.append(m * mu + b)

    # Beispielhaft nur den Fit der 302.9 keV Linie plotten
    if abs(noms(m * mu + b) - 302.85) < 2:
        K_fine = np.linspace(K_fit[0], K_fit[-1], 300)
        plt.figure()
        plt.plot(energies[K_fit], y_fit, ".", label="Data")
        plt.plot(noms(m) * K_fine + noms(b), gauss(K_fine, *p), label="Gauss fit")
        plt.xlabel("Energy $E$ / keV")
        plt.ylabel("Counts per second / s$^{-1}$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(dir + "barium_gauss_fit.pdf")
        plt.close()

N_line = np.array(N_line)
E_peak = np.array(E_peak)
for E, N in zip(E_peak, N_line):
    print(f"E = {E} keV   N = {N} 1/s")

# Emissionswahrscheinlichkeiten der Ba-133 Linien (Eckert & Ziegler Datenblatt), gleiche Reihenfolge wie peaks_fit
E_lit = np.array([276.40, 302.85, 356.02, 383.85])  # keV
W = unp.uarray([7.164, 18.33, 62.05, 8.94], [0.022, 0.06, 0.19, 0.03]) / 100

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

# Tabelle: Literaturwerte der Ba-133 Linien
tab2tex([[f"{E:.2f}" for E in E_lit], W * 100], [r"$E_\gamma$ / keV", r"$W$ / \%"], dir_tab + "ba_lit.tex",
        cpt=r"Energies $E_\gamma$ and emission probabilities $W$ of the used \ce{^{133}Ba} gamma lines~\cite{ba133}.",
        lbl="ba_lit", dgts=[0, 2])

# Tabelle: gefittete Energie, Linieninhalt, Q und Aktivität pro Linie
tab2tex([E_peak, N_line, Q, A], [r"$E$ / keV", r"$N$ / s$^{-1}$", r"$Q$", r"$A$ / Bq"], dir_tab + "ba_activity.tex",
        cpt=r"Fitted energies $E$, line contents $N$, full energy detection probabilities $Q$ and activities $A$ of the \ce{^{133}Ba} lines.",
        lbl="ba_activity", dgts=[2, 2, 2, 2])

# Diskussion: gefittete Energien im Vergleich zur Literatur
dev_ba = np.abs(noms(E_peak) - E_lit) / E_lit * 100
tab2tex([[f"{E:.2f}" for E in E_lit], E_peak, [f"{d:.3f}" for d in dev_ba]],
        [r"$E_\text{lit}$ / keV", r"$E_\text{exp}$ / keV", r"$\Delta E$ / \%"], dir_tab + "disc_ba.tex",
        cpt=r"Literature energies $E_\text{lit}$~\cite{ba133} of the \ce{^{133}Ba} lines compared with the measured energies $E_\text{exp}$.",
        lbl="disc_ba", dgts=[0, 2, 0])
