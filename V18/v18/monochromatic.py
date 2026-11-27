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
from calibration import m, b

underground = get_underground()
dir = "content/plots/"
dir_tab = "content/tables/"

# Energiekalibrierung E = m * K + b kommt aus calibration.py

# Spektrum laden: Zeile 13 bis 8204 der .Spe-Datei sind die 8192 Kanäle
counts = np.genfromtxt("raw/V18_Max_and_Julian/Caesium.Spe", skip_header=12, max_rows=8192)
counts = counts / ((1*60*60) + (18*60) + 57)  # Messzeit Caesium: 1h18min57s
counts = counts - underground

channels = np.arange(len(counts))
energies = noms(m) * channels + noms(b)  # Channels in Energie
# Kanten der Kanäle für plt.stairs (ein Wert mehr als counts)
edges = noms(m) * (np.arange(len(counts) + 1) - 0.5) + noms(b)

# Für den Fit der Photolinie (Energie, Halbwertsbreite, Zehntelwertsbreite, Linieninhalt)
def gauss(K, h, mu, sigma, c):
    return h * np.exp(-(K - mu)**2 / (2 * sigma**2)) + c

#  Compton Fit 
 
def compton(T, A_c):
    s = T / E_gamma_lit
    return A_c * (2+ s**2 / (eps**2 * (1-s)**2)+s / (1-s)*(s-2/eps))


# Theoriewerte für Compton-Kante und Rückstreupeak
E_gamma_lit = 661.657  # keV, Cs-137
m_e = const.physical_constants["electron mass energy equivalent in MeV"][0] * 1e3  # keV
eps = E_gamma_lit / m_e
E_compton_theo = E_gamma_lit * 2 * eps / (1 + 2 * eps)
E_back_theo = E_gamma_lit / (1 + 2 * eps)
print(f"Theorie: Compton-Kante = {E_compton_theo:.1f} keV, Rückstreupeak = {E_back_theo:.1f} keV")

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
print(f"Halbwertsbreite:    FWHM = {FWHM_gauss} keV")
print(f"Zehntelwertsbreite: FWTM = {FWTM_gauss} keV")

# Linieninhalt = Fläche unter dem Gauß (in Kanälen, daher Rate in 1/s)
N_photo = np.sqrt(2 * np.pi) * h * sigma
print(f"Linieninhalt Photolinie: N = {N_photo} 1/s")

K_fine = np.linspace(K_fit[0], K_fit[-1], 500)
plt.figure()
# Halbwerts- und Zehntelwertsbreite als Strecken auf der jeweiligen Höhe einzeichnen
plt.plot(noms(m) * K_fit + noms(b), y_fit, ".", label="Data")
plt.plot(noms(m) * K_fine + noms(b), gauss(K_fine, *p), label="Gauss fit")
plt.hlines(noms(c) + noms(h) / 2, noms(E_photo) - noms(FWHM_gauss) / 2, noms(E_photo) + noms(FWHM_gauss) / 2,
           color="tab:green", label="Full width at half maximum")
plt.hlines(noms(c) + noms(h) / 10, noms(E_photo) - noms(FWTM_gauss) / 2, noms(E_photo) + noms(FWTM_gauss) / 2,
           color="tab:red", label="Full width at tenth maximum")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Counts per second / s$^{-1}$")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "photopeak_Cs.pdf")
plt.close()

# Rückstreupeak: Gauß + Konstante um die erwartete Energie fitten
def fit_peak(E_guess, window):
    K_guess = int((E_guess - noms(b)) / noms(m))
    K_fit = channels[K_guess - window:K_guess + window]
    y_fit = counts[K_guess - window:K_guess + window]
    K_start = K_fit[np.argmax(y_fit)]
    p0 = [np.max(y_fit) - np.median(y_fit), K_start, window / 4, np.median(y_fit)]
    p, cov_p = curve_fit(gauss, K_fit, y_fit, p0=p0)
    h, mu, sigma, c = correlated_values(p, cov_p)
    return m * mu + b, p, K_fit, y_fit

E_back, p_back, K_back, y_back = fit_peak(E_back_theo, 200)
print(f"Rückstreupeak: gemessen = {E_back} keV, Theorie = {E_back_theo:.1f} keV")

K_fine = np.linspace(K_back[0], K_back[-1], 500)
plt.figure()
plt.plot(noms(m) * K_back + noms(b), y_back, ".", markersize=3, label="Data")
plt.plot(noms(m) * K_fine + noms(b), gauss(K_fine, *p_back), label="Gauss fit")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Counts per second / s$^{-1}$")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "backscatter_Cs.pdf")
plt.close()

# Übersichtsplot mit allen Strukturen
plt.figure()
plt.stairs(counts, edges, fill=True, label="Cs-137")
# plt.axvline(32, color="red", ls="--", label="Ba X-ray (expected)")
plt.axvline(noms(E_back), color="green", ls="--", label="Backscatter peak")
plt.axvline(E_compton_theo, color="orange", ls="--", label="Compton edge")
plt.axvline(noms(E_photo), color="pink", ls="--", label="Photopeak")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Counts per second / s$^{-1}$")
plt.yscale("log")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "monochromatic_Cs.pdf")
plt.close()


fit_mask = (energies > 350 ) & (energies < 460)
params_c, cov_c = curve_fit(compton, energies[fit_mask], counts[fit_mask], p0=[0.005])
A_c = ufloat(params_c[0], np.sqrt(cov_c[0,0]))
print(f"Compton-Fit: A_c= {A_c}")

# Integral der Fit funktion 
int_mask = (energies > 0 ) & (energies < E_compton_theo)
N_compton = A_c * np.sum(compton(energies[int_mask],1))
print(f"Linieninhalt Compton-Kontinuum: N = {N_compton} 1/s")
print(f"Verhältnis N_photo / N_compton = {N_photo / N_compton}")


E_plot = np.linspace(1, E_compton_theo - 15, 1000)

plt.figure()
plt.plot(energies[:4500], counts[:4500], ".", markersize=1, alpha=0.4, label="Data")
plt.plot(energies[fit_mask], counts[fit_mask], ".", markersize=1, color="orange", label="Fit range")
plt.plot(E_plot, compton(E_plot, noms(A_c)), color="red", label="Compton fit")
plt.xlabel("Energy $E$ / keV")
plt.ylabel("Counts per second / s$^{-1}$")
plt.legend()
plt.tight_layout()
plt.savefig(dir + "compton_Cs.pdf")
plt.close()


l = 3.9  # cm, Länge des Ge-Kristalls

# Extinktionskoeffizienten von Germanium bei 662 keV (aus extinct.png / Literatur ablesen!)  Abgelesen aus cool milimeterpapier plot
mu_photo = 0.008    # 1/cm
mu_compton = 0.37  # 1/cm

P_photo = 1 - np.exp(-mu_photo * l)
P_compton = 1 - np.exp(-mu_compton * l)
print(f"Absorptionswahrscheinlichkeit Photo: {P_photo:.4f}, Compton: {P_compton:.4f}")
print(f"Verhältnis P_photo / P_compton = {P_photo / P_compton:.4f}")