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
    return theta_v_N2 / (np.log(theta_v_N2 * R_N2 / ev + 1))


e_N2 = 3416.647339
Tv_N2 = curve_fit(e_N2)
print("Tv_N2 = ", Tv_N2)




Tv = np.linspace(300.0, 20000.0, 500)   # vibrational temperature range [K]
ev = ev_N2(Tv)
d = curve_fit(ev)

plt.figure()
plt.plot(ev, Tv)
plt.plot(ev, d)

plt.legend()
plt.ylabel(r"Vibrational Temperature $T_v$ [K]")
plt.xlabel(r"Vibrational Energy $e_v$ [J/kg]")
plt.title("N$_2$ Vibrational Energy vs Vibrational Temperature")
plt.grid(True)
plt.tight_layout()
plt.show()
