from pylab import *

data = loadtxt('ngc104.txt')
b = data[:, 0]
v = data[:, 1]
d = 4450
V = v + 5 - 5 * log10(d)

grid()
xlabel('B-V')
ylabel('$M_v$')
ylim(12,-2)
plot(b-v,V,'b.')
title('47 Tuc Cluster')
show()
