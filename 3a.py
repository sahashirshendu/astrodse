from pylab import *
from scipy.optimize import curve_fit

c=3e8
h=6.626e-34
k=1.38e-23

data=loadtxt('firas.txt')
x=data[:,0]*c*100
y=data[:,1]*1e-20
er=data[:,3]*1e-23*1000

def I(f,T):
    return 2*h/c**2*f**3/(exp(h*f/(k*T))-1)

temp,_=curve_fit(I,x,y)
print('Temperature =',temp[0],'K')
errorbar(x,y,yerr=er,label='Data')
plot(x,y,'.')
plot(x,I(x,temp[0]),label='Fitted Plot')
xlabel("$\\nu$ [Hz]")
ylabel("$I(\\nu)$ [W/$m^2$/Hz/sr]")
title('CMB')
legend()
grid()
show()
