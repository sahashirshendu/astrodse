from pylab import *

data = loadtxt('pleiades.txt')
b = data[:,0]
v = data[:,1]
d = 116
V = v + 5 - 5 * log10(d)
bv = b-v
nxl = []
nyl = []
for i in range(len(bv)):
    if (bv[i]<-0.04):
        nxl.append(bv[i])
        nyl.append(v[i])
nx = array(nxl)
ny = array(nyl)

_, pl1 = subplots()
grid()
plot(b-v,v,'b.')
plot(nx,ny,'r.')
ylabel('m_v')
xlabel('B-V')
ylim(14,2)
pl2 = pl1.twinx()
ylabel('M_v')
ylim(14+5-5*log10(d),2+5-5*log10(d))
plot(b-v,V,'b.',label='Main Sequence')
plot(nx,ny+5-5*log10(d),'r.')
legend()
title('Pleiades Cluster')
show()
