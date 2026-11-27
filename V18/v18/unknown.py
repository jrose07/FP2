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
counts = np.genfromtxt("raw/V18_Max_and_Julian/was_das_question.Spe", skip_header=12, max_rows=8192)
counts = counts / ((1*60*60) + (2*60) + 15)  # Messzeit unbekannte Quelle: 1h02min15s
counts = counts - underground

channels = np.arange(len(counts))
energies = noms(m) * channels + noms(b)  # Channels in Energie
# Kanten der Kanäle für plt.stairs (ein Wert mehr als counts)
edges = noms(m) * (np.arange(len(counts) + 1) - 0.5) + noms(b)

peaks, props = find_peaks(counts[:8000], prominence=0.05, width=4)
# Unter 150 keV: Rauschen an der Nachweisgrenze und Pb-Röntgenlinien, außerdem gilt Q dort nicht
peaks = peaks[energies[peaks] > 150]
print("Peaks in keV:", ", ".join(f"{E:.2f}" for E in energies[peaks]))

plt.figure()
plt.stairs(counts, edges, fill=True, label="Unknown source")
plt.plot(energies[peaks], counts[peaks], "x", color="red", label="peaks")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Counts per second / s$^{-1}$")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "unknown.pdf")
plt.close()

def gauss(K, h, mu, sigma, c):
    return h * np.exp(-(K - mu)**2 / (2 * sigma**2)) + c

window = 40  # Kanäle links und rechts vom Peak
N_line = []
E_peak = []
fig, axs = plt.subplots(1, len(peaks), figsize=(3 * len(peaks), 3.5))
for K, ax in zip(peaks, axs):
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
    ax.set_title(f"{noms(m * mu + b):.2f} keV")
    ax.set_xlabel("Energy $E$ / keV")
axs[0].set_ylabel("Counts per second / s$^{-1}$")
axs[0].legend()
fig.tight_layout()
fig.savefig(dir + "unknown_gauss_fits.pdf")
plt.close(fig)

N_line = np.array(N_line)
E_peak = np.array(E_peak)
for E, N in zip(E_peak, N_line):
    print(f"E = {E} keV   N = {N} 1/s")

# Zuordnung der Linien (NNDC), gleiche Reihenfolge wie peaks
nuclide = ["Ra-226", "Pb-214", "Pb-214", "Pb-214", "Bi-214"]
E_lit = np.array([186.211, 241.997, 295.224, 351.932, 609.320])  # keV
W = np.array([3.64, 7.268, 18.414, 35.60, 45.49]) / 100          # Emissionswahrscheinlichkeiten

# Q(E) = q_a * E^q_b aus calibration.py
q_a, q_b = correlated_values(params_Q, cov_Q)
Q = np.array([q_a * E**q_b for E in E_peak])

# Raumwinkel: Abstand wie bei Barium angenommen (bitte prüfen!)
a = 7.9 - 4.665 + 1.5  # Abstand Quelle - Detektor in cm
r = 4.5 / 2            # Radius des Ge-Kristalls in cm
Om = 0.5 * (1 - (a / (np.sqrt(a**2 + r**2))))  # Omega / 4pi
print("Omega/4pi =", Om)

# Aktivität A = N / (Q * W * Omega/4pi) für jede Linie
A = N_line / (Q * W * Om)
for nuc, E, q, act in zip(nuclide, E_lit, Q, A):
    print(f"{nuc}: E = {E:7.2f} keV   Q = {q}   A = {act} Bq")

# Aktivität pro Nuklid: Mittelwert über die Linien (Pb-214 hat drei)
for nuc in dict.fromkeys(nuclide):
    A_nuc = noms(A[np.array(nuclide) == nuc])
    if len(A_nuc) > 1:
        print(f"{nuc}: A = {np.mean(A_nuc):.0f} +/- {np.std(A_nuc, ddof=1) / np.sqrt(len(A_nuc)):.0f} Bq")
    else:
        print(f"{nuc}: A = {A[np.array(nuclide) == nuc][0]} Bq")

# Tabelle: gefittete Linien und Zuordnung zu den Nukliden
nuc_tex = [r"\ce{^{226}Ra}", r"\ce{^{214}Pb}", r"\ce{^{214}Pb}", r"\ce{^{214}Pb}", r"\ce{^{214}Bi}"]
tab2tex([E_peak, N_line, nuc_tex, [f"{E:.2f}" for E in E_lit], [f"{w:.2f}" for w in W * 100]],
        [r"$E$ / keV", r"$N$ / s$^{-1}$", "Nuclide", r"$E_\text{lit}$ / keV", r"$W$ / \%"], dir_tab + "unknown_lines.tex",
        cpt=r"Fitted energies $E$ and line contents $N$ of the unknown source with the assigned nuclides and their literature energies $E_\text{lit}$ and emission probabilities $W$~\cite{nndc}.",
        lbl="unknown_lines", dgts=[2, 2, 0, 0, 0])

# Tabelle: Q und Aktivität pro Linie
tab2tex([nuc_tex, [f"{E:.2f}" for E in E_lit], Q, A], ["Nuclide", r"$E_\text{lit}$ / keV", r"$Q$", r"$A$ / Bq"],
        dir_tab + "unknown_activity.tex",
        cpt=r"Full energy detection probabilities $Q$ and activities $A$ of the lines of the unknown source.",
        lbl="unknown_activity", dgts=[0, 0, 2, 2])

# Diskussion: gefittete Energien im Vergleich zur Literatur
dev_unk = np.abs(noms(E_peak) - E_lit) / E_lit * 100
tab2tex([nuc_tex, [f"{E:.2f}" for E in E_lit], E_peak, [f"{d:.3f}" for d in dev_unk]],
        ["Nuclide", r"$E_\text{lit}$ / keV", r"$E_\text{exp}$ / keV", r"$\Delta E$ / \%"], dir_tab + "disc_unknown.tex",
        cpt=r"Literature energies $E_\text{lit}$~\cite{nndc} of the assigned nuclides compared with the measured energies $E_\text{exp}$ of the unknown source.",
        lbl="disc_unknown", dgts=[0, 0, 2, 0])
