import numpy as np
import matplotlib.pyplot as plt

H = 70 / 3.08e13
wm = 1.0
a0 = 1e-10 # 0.00
a1 = 1.00
n = 100
t0 = 0


def f(x, y): # dt/da; x == a, y == t
    return x**0.5 / (H * wm**0.5)


def g(x): # t(a)
    return 2 * x**1.5 / (3 * H * wm**0.5)


def rk4(f, a, b, n, ic):
    h = (b - a) / (n - 1)
    x = np.arange(a, b + h, h)
    y = np.zeros(n)
    y[0] = ic
    for i in range(n - 1):
        y[i + 1] = y[i] + h * f(x[i] + 0.5 * h, y[i] + 0.5 * h * f(x[i], y[i]))
    return [x, y]

[a, t] = rk4(f, a0, a1, n, t0)

tr = g(a)

plt.plot(t, a, 'b-', label="Numerical")
plt.plot(tr, a, 'r.', label="Analytical")

plt.xlabel("t [seconds]")
plt.ylabel("a(t)")
plt.legend()
plt.show()
