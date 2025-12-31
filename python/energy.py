import numpy as np
import math
import matplotlib.pyplot as plt

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


Tv1 = np.linspace(300.0, 800.0, 10000)   # vibrational temperature range [K]
Tv2 = np.linspace(1950.0, 6000.0, 10000)
Tv3 = np.linspace(5950.0, 50000.0, 10000)

ev_N2_1 = ev(R_N2, Tv1, theta_v_N2)
ev_N2_2 = ev(R_N2, Tv2, theta_v_N2)
ev_N2_3 = ev(R_N2, Tv3, theta_v_N2)

ev_O2_1 = ev(R_O2, Tv1, theta_v_O2)
ev_O2_2 = ev(R_O2, Tv2, theta_v_O2)
ev_O2_3 = ev(R_O2, Tv3, theta_v_O2)

# ev_O2 = ev(R_O2, Tv, theta_v_O2)
# ev_NO = ev(R_NO, Tv, theta_v_NO)

deg = 6
N2_coeffs_1 = np.polyfit(Tv1, ev_N2_1, deg)
N2_fit_1 = np.polyval(N2_coeffs_1, Tv1)

N2_coeffs_2 = np.polyfit(Tv2, ev_N2_2, deg)
N2_fit_2 = np.polyval(N2_coeffs_2, Tv2)

N2_coeffs_3 = np.polyfit(Tv3, ev_N2_3, deg)
N2_fit_3 = np.polyval(N2_coeffs_3, Tv3)

O2_ratio_1 = ev_O2_1 / N2_fit_1

A = np.column_stack([
    np.ones_like(Tv1),
    Tv1,
    Tv1**2,
    Tv1**3,
    Tv1**4,
    np.log(Tv1),
    1.0 / np.log(Tv1)
])

coeffs = np.linalg.lstsq(A, O2_ratio_1, rcond=None)[0]
a0, a1, a2, a3, a4, b1, b2 = coeffs

O2_fit_1 = a0 + a1*Tv1 + a2*Tv1**2 + a3*Tv1**3 + a4*Tv1**4 + b1*np.log(Tv1) + b2 * 1.0 / np.log(Tv1)


plt.figure()
plt.plot(ev_N2_1, Tv1, label="N2")
plt.plot(N2_fit_1, Tv1, label="fit1", linestyle="--")
plt.plot(ev_N2_2, Tv2, label="N2")
plt.plot(N2_fit_2, Tv2, label="fit2", linestyle="--")
plt.plot(ev_N2_3, Tv3, label="N2")
plt.plot(N2_fit_3, Tv3, label="fit2", linestyle="--")


# plt.xlim(100.0, 2000.0)
plt.legend()
plt.ylabel(r"Vibrational Temperature $T_v$ [K]")
plt.xlabel(r"Vibrational Energy $e_v$ [J/kg]")
plt.title("Vibrational Energy vs Vibrational Temperature")
plt.grid(True)
plt.tight_layout()
# plt.show()



plt.figure()
plt.plot(O2_ratio_1, Tv1, label="O2/N2", linestyle="-")
# plt.plot(ev_O2_1, Tv1, label="O2")
plt.plot(O2_fit_1, Tv1, label="fit1", linestyle="--")
# plt.plot(ev_O2_2, Tv2, label="O2")
# plt.plot(O2_fit_2, Tv2, label="fit2", linestyle="--")
# plt.plot(ev_O2_3, Tv3, label="O2")
# plt.plot(O2_fit_3, Tv3, label="fit2", linestyle="--")

plt.legend()
plt.ylabel(r"Vibrational Temperature $T_v$ [K]")
plt.xlabel(r"Ratio")
plt.title("Vibrational Energy vs Vibrational Temperature")
plt.grid(True)
plt.tight_layout()
plt.show()

