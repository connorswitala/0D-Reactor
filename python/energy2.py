import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

# --- constants ---
R_u = 8314.462618          # J/mol-K

M_N2 = 28.0134         # kg/mol
R_N2 = R_u / M_N2          # J/kg-K
theta_v_N2 = 3390.0        # K

M_NO = 30.006100
R_NO = R_u / M_NO
theta_v_NO = 2739.7

M_O2 = 30.006100
R_O2 = R_u / M_O2
theta_v_O2 = 2273.54

def ev(R, T, theta_v):
    return R * theta_v / (np.exp(theta_v / T) - 1.0)


Tv1 = np.linspace(300.0, 2050.0, 10000)   # vibrational temperature range [K]
Tv2 = np.linspace(1950.0, 6000.0, 10000)
Tv3 = np.linspace(5950.0, 50000.0, 10000)

ev_N2_1 = ev(R_N2, Tv1, theta_v_N2)
ev_N2_2 = ev(R_N2, Tv2, theta_v_N2)
ev_N2_3 = ev(R_N2, Tv3, theta_v_N2)

ev_O2_1 = ev(R_O2, Tv1, theta_v_O2)
ev_O2_2 = ev(R_O2, Tv2, theta_v_O2)
ev_O2_3 = ev(R_O2, Tv3, theta_v_O2)

# Build inverse mapping: Tv(ev)
pchip = PchipInterpolator(ev_N2_1, Tv1)

# Query using ev range, not T range
x_ev = np.linspace(ev_N2_1.min(), ev_N2_1.max(), 10000)
y_Tv = pchip(x_ev)

plt.plot(ev_N2_1, Tv1, label="data (Tv vs ev)")
plt.plot(x_ev, y_Tv, "--", label="PCHIP inverse")
plt.xlabel("e_v [J/kg]")
plt.ylabel("T_v [K]")
plt.grid(True)
plt.legend()
plt.show()
