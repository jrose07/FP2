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

#Special
import io
from contextlib import redirect_stdout

with redirect_stdout(io.StringIO()):
	import stopping

dir = "content/plots/"
dir_tab = "content/tables/"

# Constants 
E_0 = 5.638*const.e*1e6 #J
m_alpha = 6.6446573450e-27 #kg
Z_air = 14.46
Z_au = 79
z = 2
rho_au = 1.93200E+01 *1e-3 / (1e-2)**3
m_atom = 196.966570 / const.N_A *1e-3 #u -> kg?

beta_0 = np.sqrt(1 - (m_alpha*const.c**2/(m_alpha*const.c**2+E_0))**2)
# print(beta_0)
#Range of alpha particles in air:
density = 1.20479E-03 #g/cm^3
range = 5.243E-03 #g/cm^2
L = range/density *1e-2 #m

"""Conversion between Voltage and Energy:"""
# STP_AIR = ...
I_air = 16 * const.e* Z_air**(0.9) #eV
I_air = 10 * const.e * Z_air
dndp = const.N_A * Z_air / ( const.R * (20+273.15)) # dn_e/dp
STP_AIR = -4*np.pi*dndp * z**2 / (const.m_e*const.c**2*beta_0**2) * (const.e**2/(4*np.pi*const.epsilon_0))**2 * ( np.log(2*const.m_e*const.c**2*beta_0**2/(I_air*(1-beta_0**2))) - beta_0**2  )

#Slopes dU/dp without foil are:
m_air = stopping.m_air #V/mbar
print(f"m_air: {m_air}")

#With that conversion factor c is:
c = m_air /L / STP_AIR #V/mbar / m /J * m * Pa = V/J * Pa/mbar
c = c * 1e-5*1e3 #Pa->mbar
#Reverse conversion (V->MeV), also E = c*U
c = 1/c /const.e /1e6 #MeV/V
print(f"1V sind {c} MeV")
c = c * const.e * 1e6 #J/V (SI units better)


"""Now calculate the thickness of the gold foil"""
dV = stopping.b_air - stopping.b_gold # Difference of energy at p = 0
dE = c * dV

#Stopping Gold:
#NIST:
STP_GOLD = 2.207e2 * const.e * 1e6 #J cm^2 / g
density_gold = 1.93200E+01 #g/cm^3
STP_GOLD *= density_gold / 1e-2 #J/m
print(STP_GOLD)
#Not NIST:
n_e = rho_au * Z_au / m_atom
I_au = 10 * const.e * Z_au
# print(n_e)
STP_GOLD = 4*np.pi*n_e * z**2 / (const.m_e*const.c**2*beta_0**2) * (const.e**2/(4*np.pi*const.epsilon_0))**2 * ( np.log(2*const.m_e*const.c**2*beta_0**2/(I_au*(1-beta_0**2))) - beta_0**2  )
print(STP_GOLD)

dx_au = dE/STP_GOLD *1e6 #µm
print(f"Die errechnete Dicke der Gold-Folie beträgt: {dx_au:.g} µm")