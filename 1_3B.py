import numpy as np
import matplotlib.pyplot as plt

H = 70 / 3.0857e19
wm = 1.0
a0 = 0.0
a1 = 1.0
n = 100

# x -> a, y -> t
def f(x, y):
    return x**0.5 / (H * wm**0.5)


def g(x):
    return 2 * x**1.5 / (3 * H * wm**0.5)

h = (a1 - a0) / (n - 1)
a = np.linspace(a0, a1, n)
t = np.zeros(n)
t[0] = 0
for i in range(n - 1):
    t[i + 1] = t[i] + h * f(a[i] + 0.5 * h, t[i] + 0.5 * h * f(a[i], t[i]))

tr = g(a)

plt.plot(t, a, 'b-', label="Numerical")
plt.plot(tr, a, 'r.', label="Analytical")

plt.xlabel("t [seconds]")
plt.ylabel("a(t)")
plt.legend()
plt.show()
