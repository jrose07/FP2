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

theta, N, dt = np.genfromtxt('raw/big_angles.csv', delimiter=",", unpack=True, skip_header=1)
#Adjust theta?
# theta = theta - (1.45 - 0.5)
Nerr = np.sqrt(N)**(-1)/dt
A = N/dt
A_err = np.sqrt(A/(dt))
tab2tex([theta, N, dt,A], [r"$\Theta \, [^{\circ}]$", r"$N$", r"$\Delta t \, [s] $", r"$A \, [s^{-1}]$"], dir_tab + "big_angles.tex", r"Measured Counts $N$ in intervall time $\Delta t$ for different angles $\Theta$, as well as the calculated count rates $A$.", lbl="big_angles", dgts=[3,5,5,3])
N = N/dt

#Save to table

def fit_func(theta, A, B):
    return A/np.pow(np.sin(np.deg2rad((theta)/2)), 4) + B

popt, pcov = curve_fit(fit_func, theta, N)

perr = np.diag(np.sqrt(pcov))
alpha = ufloat(popt[0], perr[0])
beta = ufloat(popt[1], perr[1])
print(alpha, beta)
x_fit = np.linspace(theta.min(),theta.max(),100)
y_fit = fit_func(x_fit, *popt)


fig, ax = plt.subplots()
ax.errorbar(theta, N, yerr = A_err, color="red", fmt="x", label="Data")
ax.plot(x_fit, y_fit, "b-", label=r"$\sin(\Theta/2)^{-4}$-Fit")
ax.grid()
ax.set(
    xlabel=r"$\Theta \, [°]$",
    ylabel=r"$A \, [s^{-1}]$",
    # xticks=(-4, -2, 0, 0.5, x_max, 2, 4)
    # xtickslabel=(-4,-2, 0, 0.5, )
    # ylim=[-0.5, None]
)

ax.legend()
# plt.show()
fig.savefig(dir + "big_angles.pdf")


# Now using modified fit function:

def fit_func_mod(theta, A, B, C):
    return A/np.pow(np.sin(np.deg2rad((theta - C)/2)), 4) + B

popt, pcov = curve_fit(fit_func_mod, theta, N)

perr = np.diag(np.sqrt(pcov))
alpha_mod = ufloat(popt[0], perr[0])
beta_mod = ufloat(popt[1], perr[1])
theta_0_mod = ufloat(popt[2], perr[2])
print(alpha, beta, theta_0_mod)
x_fit = np.linspace(theta.min(),theta.max(),100)
y_fit = fit_func_mod(x_fit, *popt)


fig, ax = plt.subplots()
ax.errorbar(theta, N, yerr = A_err, color="red", fmt="x", label="Data")
ax.plot(x_fit, y_fit, "b-", label=r"$\sin((\Theta - \Theta_0)/2)^{-4}$-Fit")
ax.grid()
ax.set(
    xlabel=r"$\Theta \, [°]$",
    ylabel=r"$A \, [s^{-1}]$",
    # xticks=(-4, -2, 0, 0.5, x_max, 2, 4)
    # xtickslabel=(-4,-2, 0, 0.5, )
    # ylim=[-0.5, None]
)

ax.legend()
# plt.show()
fig.savefig(dir + "big_angles_mod.pdf")



#Fit parameters 
write('Fit Parameters for Rutherford cross section:\n')
add(f"alpha = {alpha:.2g} s^(-1), beta = {beta:.2g} s^(-1)\n")
add('Fit Parameters for modified fit:\n')
add(f"alpha = {alpha_mod:.2g} s^(-1), beta = {beta_mod:.2g} s^(-1), Theta_0 = {theta_0_mod:.2g} DEG\n\n")