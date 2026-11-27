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
from underground import get_underground

underground = get_underground()
dir = "content/plots/"
dir_tab = "content/tables/"
t = 824428800
A_0 = ufloat(4130,60)
T_2 = 13.537*365.25*24*60*60
a =  7.9-4.675 +1.5 # Abstand Quelle - Detektor in cm
r = 4.5 / 2      # Radius des Ge-Kristalls in cm (Durchmesser 45 mm)

# Spektrum laden: Zeile 13 bis 8204 der .Spe-Datei sind die 8192 Kanäle 
counts = np.genfromtxt("raw/V18_Max_and_Julian/EU_Callibrate.Spe", skip_header=12, max_rows=8192)
counts = counts / ((1*60*60) + (35))
counts = counts - underground


channels = np.arange(len(counts))

# width filtert einzelne Rausch-Kanäle raus, Kanäle ab 8000 sind Overflow
peaks, props = find_peaks(counts[:8000], prominence=0.005, width=4)

# Eu-152 Linien aus der Tabelle, die im Messbereich liegen (> ~950 keV ist im Overflow)
E_lit = np.array([121.7817, 244.6975, 344.2785, 411.1163, 443.965, 778.9040, 867.378])
# Emissionswahrscheinlichkeiten in % aus dem Eckert & Ziegler Datenblatt (mit Unsicherheit)
I_lit = unp.uarray([28.58, 7.580, 26.5, 2.234, 3.148, 12.94, 4.245],
                   [0.09, 0.030, 0.6, 0.025, 0.020, 0.15, 0.021])

# Grobe Zuordnung über die zwei stärksten Peaks (121.78 keV und 344.28 keV)
strongest = np.sort(peaks[np.argsort(props["prominences"])[-2:]])
m_grob = (E_lit[2] - E_lit[0]) / (strongest[1] - strongest[0])
b_grob = E_lit[0] - m_grob * strongest[0]

# Zu jeder Literaturenergie den nächsten gefundenen Peak nehmen
K_peak = np.array([peaks[np.argmin(np.abs(peaks - (E - b_grob) / m_grob))] for E in E_lit])
print("Zugeordnete Kanäle:", K_peak)

# Tabelle: Literaturwerte der Eu-152 Gammalinien und die zugeordneten Kanäle
E_tab = ["121.78", "244.70", "344.28", "411.12", "443.96", "778.90", "867.37"]
tab2tex([E_tab, I_lit, [str(K) for K in K_peak]],
        [r"$E_\gamma$ / keV", r"$W$ / \%", r"Channel $K$"], dir_tab + "eu_lit.tex",
        cpt=r"Energies $E_\gamma$ and emission probabilities $W$ of the used \ce{^{152}Eu} gamma lines~\cite{eu152} and the assigned channels $K$ of the measured peaks.",
        lbl="Euspec", dgts=[0, 2, 0])

A = A_0*(np.exp(-(np.log(2)/T_2)*t))
print("A=",A)

Om = 0.5* (1-(a/(np.sqrt(a**2+r**2))))  # Omega / 4pi
print("Omega/4pi =", Om)



plt.figure()
plt.stairs(counts, np.arange(len(counts) + 1), fill=True, label="Eu-152")
plt.plot(K_peak, counts[K_peak], "x", color="red", label="peaks")
plt.xlabel("Channel")
plt.ylabel("Counts per second / s$^{-1}$")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "eu_calibration.pdf")
plt.close()

# Lineare Ausgleichsrechnung E = m * K + b
def linear(K, m, b):
    return m * K + b

params, cov = curve_fit(linear, K_peak, E_lit)
errors = np.sqrt(np.diag(cov))
m = ufloat(params[0], errors[0])
b = ufloat(params[1], errors[1])
print(f"m = {m} keV/Kanal")
print(f"b = {b} keV")

K_plot = np.linspace(0, 8000, 2)
plt.figure()
plt.plot(K_peak, E_lit, "x", label="Eu-152 peaks")
plt.plot(K_plot, linear(K_plot, *params), label="Linear fit")
plt.xlabel("Channel $K$")
plt.ylabel("Energy $E$ / keV")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "eu_linear_fit.pdf")
plt.close()

# Linieninhalt: Gauß + Konstante (Untergrund unter dem Peak) an jeden Peak fitten
from uncertainties import correlated_values

def gauss(K, h, mu, sigma, c):
    return h * np.exp(-(K - mu)**2 / (2 * sigma**2)) + c

window = 40  # Kanäle links und rechts vom Peak
N_line = []
E_fit = []
for K, E in zip(K_peak, E_lit):
    K_fit = channels[K - window:K + window]
    y_fit = counts[K - window:K + window]
    p0 = [counts[K], K, 5, np.median(y_fit)]
    p, cov_g = curve_fit(gauss, K_fit, y_fit, p0=p0)
    h, mu, sigma, c = correlated_values(p, cov_g)
    # Fläche unter dem Gauß = Linieninhalt (hier als Rate in 1/s, weil counts schon durch t geteilt ist)
    N_line.append(np.sqrt(2 * np.pi) * h * abs(sigma))
    E_fit.append(m * mu + b)

    # Beispielhaft nur den Fit der 344.3 keV Linie plotten
    if np.isclose(E, 344.2785):
        K_fine = np.linspace(K_fit[0], K_fit[-1], 300)
        plt.figure()
        plt.plot(K_fit, y_fit, ".", label="Data")
        plt.plot(K_fine, gauss(K_fine, *p), label="Gauss fit")
        plt.xlabel("Channel $K$")
        plt.ylabel("Counts per second / s$^{-1}$")
        plt.legend()
        plt.tight_layout()
        plt.savefig(dir + "eu_gauss_fit.pdf")
        plt.close()
N_line = np.array(N_line)

# Full energy detection probability Q = N / (A * W * Omega/4pi)
Q = N_line / (A * I_lit / 100 * Om)
for E, N, q in zip(E_lit, N_line, Q):
    print(f"E = {E:8.2f} keV   N = {N} 1/s   Q = {q}")

# Tabelle: Linieninhalt und Vollenergienachweiswahrscheinlichkeit pro Linie
tab2tex([E_tab, N_line, Q], [r"$E_\gamma$ / keV", r"$N$ / s$^{-1}$", r"$Q$"], dir_tab + "eu_Q.tex",
        cpt=r"Line contents $N$ and full energy detection probabilities $Q$ of the \ce{^{152}Eu} lines.",
        lbl="eu_Q", dgts=[0, 2, 2])

# Potenzfunktion nur für E > 150 keV fitten (siehe Anleitung)
def power(E, q_a, q_b):
    return q_a * E**q_b

mask = E_lit > 150
params_Q, cov_Q = curve_fit(power, E_lit[mask], noms(Q[mask]), sigma=stds(Q[mask]), absolute_sigma=True, p0=[30, -1])
q_a, q_b = correlated_values(params_Q, cov_Q)
print(f"Q(E) = a * E^b mit a = {q_a}, b = {q_b}")

E_plot = np.linspace(150, 950, 200)
plt.figure()
plt.errorbar(E_lit[mask], noms(Q[mask]), yerr=stds(Q[mask]), fmt="x", label="Data (fit)")
plt.errorbar(E_lit[~mask], noms(Q[~mask]), yerr=stds(Q[~mask]), fmt="x", color="gray", label="Data (not used)")
plt.plot(E_plot, power(E_plot, *params_Q), label="Power fit")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Full energy detection probability $Q$")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "eu_efficiency.pdf")
plt.close()

# Diskussion: Energie der Peaks aus der Kalibrierung im Vergleich zur Literatur
E_fit = np.array(E_fit)
dev_eu = np.abs(noms(E_fit) - E_lit) / E_lit * 100
tab2tex([E_tab, E_fit, [f"{d:.3f}" for d in dev_eu]],
        [r"$E_\text{lit}$ / keV", r"$E_\text{exp}$ / keV", r"$\Delta E$ / \%"], dir_tab + "disc_eu.tex",
        cpt=r"Literature energies $E_\text{lit}$~\cite{eu152} of the \ce{^{152}Eu} lines compared with the energies $E_\text{exp}$ from the gaussian fits and the calibration.",
        lbl="disc_eu", dgts=[0, 2, 0])
