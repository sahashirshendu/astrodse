from pylab import *

H=70/3.086e19
a0=0.0
a1=1.0
n=1000

def f(a,t):
    return a**0.5/H

def g(a):
    return 2*a**1.5/(3*H)

h=(a1-a0)/(n-1)
a=linspace(a0,a1,n)
t=zeros(n)
t[0]=0
for i in range(n - 1):
    t[i+1]=t[i]+h*f(a[i]+.5*h,t[i]+.5*h*f(a[i],t[i]))

plot(t,a,'-',label='Numerical')
plot(g(a),a,'.',label='Analytical')
xlabel('t [seconds]')
ylabel('a(t)')
legend()
show()
