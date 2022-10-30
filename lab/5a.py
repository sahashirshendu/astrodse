from pylab import *

data = loadtxt('pleiades.txt')
b = data[:, 0]
v = data[:, 1]
bv = b-v
nx = []
ny = [] # not in main sequenc
for i in range(len(bv)):
    if (bv[i]<-0.04):
        nx.append(bv[i])
        ny.append(v[i])

grid()
xlabel('B-V')
ylabel('m_v')
ylim(14,2)
plot(b-v,v,'b.',label='Main Sequence')
plot(nx,ny,'r.')
title('Pleiades Cluster')
legend()
show()
