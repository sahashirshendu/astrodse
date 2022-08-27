import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({"text.usetex": True})

H = 70 / 3.0857e19
wm = 0.3
wr = 5e-5
wde = 0.7
w = -1.0
a0 = 0.0
a1 = 1.0
n = 100

# x -> a, y -> t
def f(x, y):
    return 1 / (H * np.sqrt(wm / x + wr / x ** 2 + wde / x ** (3 * w - 1)))

h = (a1 - a0) / (n - 1)
a = np.linspace(a0, a1, n)
t = np.zeros(n)
t[0] = 0
for i in range(n - 1):
    t[i + 1] = t[i] + h * f(a[i] + 0.5 * h, t[i] + 0.5 * h * f(a[i], t[i]))

plt.plot(t, a, 'b-', label="Numerical")
plt.xlabel("t [seconds]")
plt.ylabel("a(t)")
plt.legend()
plt.show()
