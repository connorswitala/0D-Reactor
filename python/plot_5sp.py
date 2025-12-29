import numpy as np
import matplotlib.pyplot as plt
import subprocess

subprocess.run(["cmake", "--build", "."], cwd="../build", check=True)
subprocess.run(["./../build/reactor.exe"], check=True)

data = np.genfromtxt("../files/0Dreactor.csv", delimiter=",", names=True)

t     = data["t"]
T_tr  = data["T_tr"]
T_v   = data["T_v"]

X_N2  = data["X_N2"]
X_N   = data["X_N"]
X_O2  = data["X_O2"]
X_O   = data["X_O"]
X_NO  = data["X_NO"]

# Mask zeros for log plots
X_NO[X_NO <= 0] = np.nan

# --------------------
# Figure 1: Temperatures
# --------------------
plt.figure()
plt.plot(t, T_tr, label="T_tr")
plt.plot(t, T_v,  label="T_v")

plt.xscale("log")
plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")

plt.minorticks_on()
plt.grid(which="major", linestyle="-",  alpha=0.7)
plt.grid(which="minor", linestyle="--", alpha=0.4)

plt.legend()
plt.tight_layout()
plt.savefig("park_temp.png", dpi=300)
plt.close()

# --------------------
# Figure 2: Species
# --------------------
plt.figure()
plt.plot(t, X_N2, label="N2")
plt.plot(t, X_N,  label="N")
plt.plot(t, X_O2, label="O2")
plt.plot(t, X_O,  label="O")
plt.plot(t, X_NO, label="NO")

plt.xscale("log")

plt.xlabel("Time [s]")
plt.ylabel("Mole fraction")

plt.minorticks_on()
plt.grid(which="major", linestyle="-",  alpha=0.7)
plt.grid(which="minor", linestyle="--", alpha=0.4)

plt.legend()
plt.tight_layout()
plt.savefig("park_species.png", dpi=300)
plt.close()
