from pylab import *
from scipy.integrate import quad

H_0=70/3.086e19
c=3e8
dh=(c/H_0)/3.086e22
th=1/H_0
wm=0.3
wr=5e-5
wl=0.7
w=-1.0
z=5.0
n=1000

h=z/(n-1)
da=zeros(n)
dl=zeros(n)

def EO(z):
    return (wm*(1+z)**3 + wr*(1+z)**4 + wl*(1+z)**(3*(1+w)))**-0.5

def ET(z):
    return (wm*(1+z)**3 + wr*(1+z)**4 + wl*(1+z)**(3*(1+w)))**-0.5 / (1+z)

zt=10000
nt=100000
ht=zt/(nt-1)

zpe=0
tsum=ET(0) + ET(zt)
for i in range (1,nt-1):
    zpe = i*ht
    tsum = tsum + 2 * ET(zpe)
t = th*(ht/2*tsum)/(pi*1e16)

zre=0
dpe=0
zr=[]
dp=[]
dsum=0
for i in range (1,n+1):
    zr.append(zre)
    dp.append(dh * quad(EO,0,zre)[0])
    zre=i*h

for i in range(n):
    da[i]=dp[i]/(1+zr[i])
    dl[i]=dp[i]*(1+zr[i])

print("Age of the universe =",t,"billion years")
plot(zr,dp,label="Distance")
plot(zr,da,label="Angular Diameter Distance")
plot(zr,dl,label="Luminosity Distance")
xlabel("Redshift (z)")
ylabel("Distance (Megaparsec)")
legend()
grid()
show()
