import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt("pleiades.txt")
b = data[:, 0]
v = data[:, 1]
d = 116
V = v + 5 - 5 * np.log10(d)

fig, ax1 = plt.subplots()

ax1.grid()
ax1.set_xlabel("$B-V$")
ax1.set_ylabel("$m_V$")
ax1.invert_yaxis()
ax1.set_xlim(-0.5,2)
# ax1.set_ylim(max(MV), min(MV))
ax1.plot(b - v, v, "b.")
ax1.set_xticks([-0.35, -0.31, -0.16, 0.00, 0.13, 0.27, 0.42, 0.58, 0.70, 0.89, 1.18, 1.45, 1.63])
ax1.set_xticklabels(["O5", "B0", "B5", "A0", "A5", "F0", "F5", "G0", "G5", "K0", "K5", "M0", "M5"])
 
ax2 = ax1.twinx()
ax2.set_ylabel("$M_V$")
ax2.invert_yaxis()
# ax2.set_ylim(max(MV), min(MV))
ax2.plot(b - v, V, "b.")
ax2.scatter(0.70, 4.83, marker="*", label="Sun")
ax2.legend()

ax3 = ax1.twiny()
ax3.plot(b - v, v, "b.")

plt.title("Pleiades")
# plt.grid()
# plt.gca().invert_yaxis()
plt.show()
