import numpy as np
import math
import matplotlib.pyplot as plt

# --- constants ---
R_u = 8.314462618          # J/mol-K
M_N2 = 28.0134e-3          # kg/mol
R_N2 = R_u / M_N2          # J/kg-K
theta_v_N2 = 3390.0        # K

def ev_N2(Tv):
    return R_N2 * theta_v_N2 / (np.exp(theta_v_N2 / Tv) - 1.0)

def curve_fit(ev):
    return 1 / (R_N2 - 2.2) * ev + 1500

def curve_fit2(ev):
    return theta_v_N2 / (np.log(theta_v_N2 * R_N2 / ev + 1))


print(1 / R_N2)

Tv = np.linspace(300.0, 20000.0, 500)   # vibrational temperature range [K]
ev = ev_N2(Tv)
c = curve_fit(ev)
d = curve_fit2(ev)

plt.figure()
plt.plot(ev, Tv)
plt.plot(ev, c)
plt.plot(ev, d)

plt.legend()
plt.ylabel(r"Vibrational Temperature $T_v$ [K]")
plt.xlabel(r"Vibrational Energy $e_v$ [J/kg]")
plt.title("N$_2$ Vibrational Energy vs Vibrational Temperature")
plt.grid(True)
plt.tight_layout()
plt.show()
