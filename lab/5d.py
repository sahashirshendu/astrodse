from pylab import *

data = loadtxt("pleiades.txt")
b = data[:, 0]
v = data[:, 1]
d = 116
V = v + 5 - 5 * log10(d)

grid()
xlabel("B-V")
ylabel("$M_v$")
ylim(8,-4)
plot(b-v,V,"b.")
title("Pleiades Cluster")
show()
