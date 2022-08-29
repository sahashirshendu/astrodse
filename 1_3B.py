import numpy as np
import matplotlib.pyplot as plt

H = 70 / 3.086e19
a0 = 0.0
a1 = 1.0
n = 100

def f(a, t):
    return a**0.5 / H

def g(a):
    return 2 * a**1.5 / (3 * H)

h = (a1 - a0) / (n - 1)
a = np.linspace(a0, a1, n)
t = np.zeros(n)
t[0] = 0
for i in range(n - 1):
    t[i + 1] = t[i] + h * f(a[i] + 0.5 * h, t[i] + 0.5 * h * f(a[i], t[i]))

plt.plot(t,a,'-',label='Numerical')
plt.plot(g(a),a,'.',label='Analytical')
plt.xlabel('t [seconds]')
plt.ylabel('a(t)')
plt.legend()
plt.show()
