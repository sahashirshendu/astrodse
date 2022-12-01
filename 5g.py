from pylab import *

data = loadtxt('ngc104.txt')
b = data[:,0]
v = data[:,1]
d = 4450
V = v + 5 - 5 * log10(d)
wd = loadtxt('wd.txt')
wx = wd[:,0]
wy = wd[:,1]
bv = b-v
rxl = []
ryl = []
for i in range(len(bv)):
    if (bv[i]>0.55) and (v[i]<17):
        rxl.append(bv[i])
        ryl.append(v[i])
rx = np.array(rxl)
ry = np.array(ryl)

_, pl1 = subplots()
plot(b-v,v,'b.')
plot(wx,wy,'g.')
plot(rx,ry,'r.')
xlabel('B-V')
ylabel('m_v')
ylim(25,10)
pl2 = pl1.twinx()
ylabel('M_v')
ylim(25+5-5*log10(d),10+5-5*log10(d))
plot(b-v,V,'b.',label='Main Sequence')
plot(wx,wy+5-5*log10(d),'g.',label='White Dwarfs')
plot(rx,ry+5-5*log10(d),'r.',label='Red Giants')
title('47 Tuc Cluster')
legend()
show()
