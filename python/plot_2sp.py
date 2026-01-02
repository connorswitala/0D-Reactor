import numpy as np
import matplotlib.pyplot as plt
import subprocess

subprocess.run(["cmake", "--build", "."], cwd="../build", check=True)
subprocess.run(["./../build/N2reactor.exe"], check=True)

data = np.genfromtxt("../files/0Dreactor.csv", delimiter=",", names=True)

t     = data["t"]
T_tr  = data["T_tr"]
T_v   = data["T_v"]

X_N2  = data["X_N2"]
X_N   = data["X_N"]


# --------------------
# Figure 1: Temperatures
# --------------------
plt.figure()
plt.plot(t, T_tr, label="T_tr")
plt.plot(t, T_v,  label="T_v")

plt.xscale("log")
plt.xlim(1e-10, 1e-6)

plt.xlabel("Time [s]")
plt.ylabel("Temperature [K]")

plt.minorticks_on()
plt.grid(which="major", linestyle="-",  alpha=0.7)
plt.grid(which="minor", linestyle="--", alpha=0.4)

plt.legend()
plt.tight_layout()
plt.savefig("N2_temperature_vs_time2.png", dpi=300)
plt.close()

# --------------------
# Figure 2: Species
# --------------------
plt.figure()
plt.plot(t, X_N2, label="N2")
plt.plot(t, X_N,  label="N")

plt.xscale("log")
plt.xlim(1e-10, 1e-6)


plt.xlabel("Time [s]")
plt.ylabel("Mole fraction")

plt.minorticks_on()
plt.grid(which="major", linestyle="-",  alpha=0.7)
plt.grid(which="minor", linestyle="--", alpha=0.4)

plt.legend()
plt.tight_layout()
plt.savefig("N2_species_vs_time2.png", dpi=300)
plt.close()
