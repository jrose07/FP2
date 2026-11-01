import uncertainties.unumpy as unp
import numpy as np
from scipy.optimize import curve_fit
from addons import write, add, latex_float, tab_to_latex as tab2tex
from scipy.stats import linregress
from uncertainties import ufloat
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
import scipy.constants as const
import matplotlib.pyplot as plt

dir = "content/plots/"
dir_tab = "content/tables/"

theta, N = np.genfromtxt('raw/small_angles.csv', delimiter=",", unpack=True, skip_header=1)
N = N/20
# print(theta, N)

# def fit_func(theta, A, B):
#     return A/np.pow(np.sin(np.deg2rad(theta)), 4) + B

def fit_func(theta, a, x0, sigma):
# def gaussian(x, a, x0, sigma):
    return a * np.exp(-(theta - x0)**2 / (2 * sigma**2))

popt, pcov = curve_fit(fit_func, theta, N)

perr = np.diag(np.sqrt(pcov))
print(perr)
x_fit = np.linspace(theta.min(),theta.max(),100)
y_fit = fit_func(x_fit, *popt)

x_max = x_fit[np.argmax(y_fit)]
print(x_max)
fig, ax = plt.subplots()
ax.errorbar(theta, N, yerr = 20/np.sqrt(N*20), color="red", fmt="x", label="Data")
ax.plot(x_fit, y_fit, "b-", label="Gaussian Fit")
# ax.vlines(x_max, ymin=0*0.9*y_fit.max(), ymax=1.1*y_fit.max(), color="orange", linestyles="dashed")
selected = theta == 0.5
ax.axvline(x_max, color="orange", linestyle="dashed", label=r"optimal $\Theta_0$")
ax.axvline(0.5, color="green", linestyle="dashed", label=r"selected $\Theta_0$")
ax.grid()
ax.set(
    xlabel=r"$\Theta \, [°]$",
    ylabel=r"$N \, [s^{-1}]$",
    xticks=(-4, -2, 0, 0.5, x_max, 2, 4)
    # xtickslabel=(-4,-2, 0, 0.5, )
    # ylim=[-0.5, None]
)
ax.set_xticklabels([-4, -2, 0, 0.5, f"{x_max:.2f}", 2,4])

ax.legend()
# plt.show()
fig.savefig(dir + "small_angles.pdf")