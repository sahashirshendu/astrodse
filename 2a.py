from pylab import *

H_0=70/3.086e19
c=3e8
dh=(c/H_0)/3.086e22
th=1/H_0
wm=0.0
wr=1./3.
z=5.0
n=1000

zp=linspace(0,z,n)
dpm=zeros(n)
dam=zeros(n)
dlm=zeros(n)
dpr=zeros(n)
dar=zeros(n)
dlr=zeros(n)

for i in range(n):
    dpm[i] = 2*dh/(3*wm+1)*(1-(1+zp[i])**(-(3*wm+1)/2))
    dpr[i] = 2*dh/(3*wr+1)*(1-(1+zp[i])**(-(3*wr+1)/2))
    dam[i] = dpm[i]/(1+zp[i])
    dar[i] = dpr[i]/(1+zp[i])
    dlm[i] = dpm[i]*(1+zp[i])
    dlr[i] = dpr[i]*(1+zp[i])

agem=(th*2/(3*(1+wm)))/(pi*1e16)
ager=(th*2/(3*(1+wr)))/(pi*1e16)

print("Age of matter only universe =",agem,"billion years")
print("Age of radiation only universe =",ager,"billion years")
plot(zp,dpm,label="Distance (Matter Only)")
plot(zp,dam,label="Angular Diameter Distance (Matter Only)")
plot(zp,dlm,label="Luminosity Distance (Matter Only)")
plot(zp,dpr,label="Distance (Radiation Only)")
plot(zp,dar,label="Angular Diameter Distance (Radiation Only)")
plot(zp,dlr,label="Luminosity Distance (Radiation Only)")
xlabel("Redshift (z)")
ylabel("Distance (Megaparsec)")
legend()
grid()
show()
