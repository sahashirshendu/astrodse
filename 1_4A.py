import numpy as np
import matplotlib.pyplot as plt

H_0 = 70/3.086e19
c = 3e8
dh = (c/H_0) / 3.086e22
th = 1/H_0

wm = 0.0
wr = 1./3.
z = 5.0
n = 1000

zp = np.linspace(0, z, n)
dpm = np.zeros(n)
dam = np.zeros(n)
dlm = np.zeros(n)
dpr = np.zeros(n)
dar = np.zeros(n)
dlr = np.zeros(n)

def dpp(z, w):
    return 2 * dh / (3 * w + 1) * (1 - (1 + z) ** (- (3 * w + 1) / 2))

for i in range(n):
    dpm[i] = dpp(zp[i], wm)
    dpr[i] = dpp(zp[i], wr)
    dam[i] = dpm[i] / (1 + zp[i])
    dar[i] = dpr[i] / (1 + zp[i])
    dlm[i] = dpm[i] * (1 + zp[i])
    dlr[i] = dpr[i] * (1 + zp[i])

agem = (th * 2 / (3*(1+wm))) / (np.pi*1e16)
ager = (th * 2 / (3*(1+wr))) / (np.pi*1e16)

print("Age of matter only universe =",agem,"billion years")
print("Age of radiation only universe =",ager,"billion years")
plt.plot(zp,dpm,label="Proper Distance (Matter Only)")
plt.plot(zp,dam,label="Angular Diameter Distance (Matter Only)")
plt.plot(zp,dlm,label="Luminosity Distance (Matter Only)")
plt.plot(zp,dpr,label="Proper Distance (Radiation Only)")
plt.plot(zp,dar,label="Angular Diameter Distance (Radiation Only)")
plt.plot(zp,dlr,label="Luminosity Distance (Radiation Only)")
plt.xlabel("Redshift (z)")
plt.ylabel("Distance (Megaparsec)")
plt.legend()
plt.grid()
plt.show()
