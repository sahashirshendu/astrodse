from pylab import *

data = loadtxt('pleiades.txt')
b = data[:, 0]
v = data[:, 1]
d = 116
V = v + 5 - 5 * log10(d)

_, pl1 = subplots()
grid()
plot(b-v,v,'b.')
ylabel('m_v')
ylim(14,2)
pl2 = pl1.twinx()
xlabel('B-V')
ylabel('M_v')
ylim(14+5-5*log10(d),2+5-5*log10(d))
plot(b-v,V,'b.')
scatter(0.656,4.74,marker='*',label='Sun')
title('Pleiades Cluster')
show()
