import uncertainties.unumpy as unp
import numpy as np
from scipy.optimize import curve_fit
from addons import write, add, latex_float, tab_to_latex as tab2tex
from scipy.stats import linregress
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
import scipy.constants as const
import matplotlib.pyplot as plt
import pandas as pd

dir = "content/plots/"
dir_tab = "content/tables/"

#Air
p, dV = np.genfromtxt("raw/stopping_air.csv", delimiter=",", skip_header=1, unpack=True)

fig, ax = plt.subplots()
ax.grid()

#Select the ones outside of Uderground
selected = np.logical_and(dV >= 5.8, p<=400)
notselect = np.logical_not(selected)
ax.errorbar(p[notselect], dV[notselect], xerr=1, yerr=0.01, color="red", fmt="x", label="Underground Data")
p_air, dV_air = p[selected], dV[selected]
ax.errorbar(p_air, dV_air, xerr=1, yerr=0.01, color="green", fmt="x", label="Selected Data")

#Linregress
result = linregress(p_air, dV_air)
m = ufloat(result.slope, result.stderr)
b = ufloat(result.intercept, result.intercept_stderr)
p_fita = np.linspace(p_air.min()-10, p.max(), 100)
dV_fita = noms(m)*p_fita + noms(b)
ax.plot(p_fita, dV_fita, "b-", label="Linear Regression")
m_air = m
b_air = b

ax.set(
    xlabel=r"$p \, [\text{mbar}]$",
    ylabel=r"$\Delta V \, [V]$",
    ylim=[0.9*dV.min(), 1.1*dV.max()]
)
ax.legend()
fig.savefig(dir + "stopping_air.pdf")

#Gold
p, dV = np.genfromtxt("raw/stopping_gold.csv", delimiter=",", skip_header=1, unpack=True)

fig, ax = plt.subplots()
ax.grid()

#Select the ones outside of Uderground
selected = np.logical_and(dV >= 5.8, p<=400)
selected = np.logical_and(dV >= 5.8, p<=400)
notselect = np.logical_not(selected)
ax.errorbar(p[notselect], dV[notselect], xerr=1, yerr=0.01, color="red", fmt="x", label="Underground Data")
p_gold, dV_gold = p[selected], dV[selected]
ax.errorbar(p_gold, dV_gold, xerr=1, yerr=0.01, color="green", fmt="x", label="Selected Data")


#Linregress
result = linregress(p_gold, dV_gold)
m = ufloat(result.slope, result.stderr)
b = ufloat(result.intercept, result.intercept_stderr)
m_gold = m
b_gold = b
p_fitg = np.linspace(p_gold.min()-20, p.max(), 100)
dV_fitg = noms(m)*p_fitg + noms(b)
ax.plot(p_fitg, dV_fitg, "b-", label="Linear Regression")

ax.set(
    xlabel=r"$p \, [\text{mbar}]$",
    ylabel=r"$\Delta V \, [V]$",
    ylim=[0.9*dV.min(), 1.1*dV.max()],
    xlim=[-50+p.min(), 1.02*p.max()]
)
ax.legend()
fig.savefig(dir + "stopping_gold.pdf")

#Together
fig, ax = plt.subplots()
ax.plot(p_air, dV_air, "rx", label=r"Data without foil")
ax.plot(p_fita, dV_fita, color="orange", linestyle="-", label=r"Linear Regression without foil")
ax.plot(p_gold, dV_gold, "gx", label="Data with foil")
ax.plot(p_fitg, dV_fitg, color="limegreen", linestyle="-", label="Linear Regression with foil")
ax.set(
    xlabel=r"$p \, [\text{mbar}]$",
    ylabel=r"$\Delta V \, [V]$",
    ylim=[-0.1*dV.min(), 1.2*dV.max()],
    xlim=[-50+p_air.min(), 2.1*p_air.max()]
)
ax.grid()
ax.legend()
# plt.show()
fig.savefig(dir + "stopping_comparison.pdf")

"""Calculate the dx of the gold foil now"""
print(m_gold, m_air)
print(b_gold, b_air)

add('Linear Regressions of Stopping:\n')
add(f"With Foil:\nm={m_gold:.4f}V/mbar, b={b_gold:.4f}mbar\n")
add(f"Without Foil:\nm={m_air:.4f}V/mbar, b={b_air:.4f}mbar\n")