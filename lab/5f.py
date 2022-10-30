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
scatter(0.656,4.74,marker='*',label='Sun')
xticks([-.35,-.31,-.16,0,.13,.27,.42,.58,.7,.89,1.18,1.45,1.63],['O5','B0','B5','A0','A5','F0','F5','G0','G5','K0','K5','M0','M5'])
title('Pleiades Cluster')
legend()
show()
