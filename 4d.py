import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

H_0 = 70/3.086e19
c = 3e8
dh = (c/H_0) / 3.086e22
th = 1/H_0

wm = 0.3
wr = 5e-5
wl = 0.7
w = -1.0
w1 = -2.0/3.0
w2 = -1.0/2.0
z = 5.0
n = 1000

h = z / (n-1)
da1 = np.zeros(n)
dl1 = np.zeros(n)
da2 = np.zeros(n)
dl2 = np.zeros(n)

def EO1(z):
    return (wm*(1+z)**3 + wr*(1+z)**4 + wl*(1+z)**(3*(1+w1)))**-0.5

def EO2(z):
    return (wm*(1+z)**3 + wr*(1+z)**4 + wl*(1+z)**(3*(1+w2)))**-0.5

def ET(z, w):
    return 1 / (1+z) * (wm*(1+z)**3 + wr*(1+z)**4 + wl*(1+z)**(3*(1+w)))**-0.5

zt = 10000
nt = 100000
ht = zt / (nt-1)

zpe = 0
tsum1 = ET(0, w1) + ET(zt, w1)
tsum2 = ET(0, w2) + ET(zt, w2)
for i in range (1,nt-1):
    zpe = i*ht
    tsum1 = tsum1 + 2 * ET(zpe, w1)
    tsum2 = tsum2 + 2 * ET(zpe, w2)
t1 = th * (ht/2 * tsum1) / (np.pi*1e16)
t2 = th * (ht/2 * tsum2) / (np.pi*1e16)

zre = 0
dpe1 = 0
dpe2 = 0
zr = []
dp1 = []
dp2 = []
dsum1 = 0
dsum2 = 0
for i in range (1,n+1):
    zr.append(zre)
    dp1.append(dh * quad(EO1, 0, zre)[0])
    dp2.append(dh * quad(EO2, 0, zre)[0])
    zre = i*h

for i in range(n):
    da1[i] = dp1[i] / (1 + zr[i])
    dl1[i] = dp1[i] * (1 + zr[i])
    da2[i] = dp2[i] / (1 + zr[i])
    dl2[i] = dp2[i] * (1 + zr[i])

print("Age of universe (w=-2/3) =",t1,"billion years")
print("Age of universe (w=-1/2) =",t2,"billion years")
plt.plot(zr,dp1,label="Proper Distance (w = -2/3)")
plt.plot(zr,da1,label="Angular Diameter Distance (w = -2/3)")
plt.plot(zr,dl1,label="Luminosity Distance (w = -2/3)")
plt.plot(zr,dp2,label="Proper Distance (w = -1/2)")
plt.plot(zr,da2,label="Angular Diameter Distance (w = -1/2)")
plt.plot(zr,dl2,label="Luminosity Distance (w = -1/2)")
plt.xlabel("Redshift (z)")
plt.ylabel("Distance (Megaparsec)")
plt.legend()
plt.grid()
plt.show()
