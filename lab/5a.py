from pylab import *

data = loadtxt("pleiades.txt")
b = data[:, 0]
v = data[:, 1]

grid()
xlabel("B-V")
ylabel("$m_v$")
ylim(14,2)
plot(b-v,v,"b.")
title("Pleiades Cluster")
show()
