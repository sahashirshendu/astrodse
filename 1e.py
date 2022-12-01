from pylab import *

H=70/3.086e19
wr=5e-5
wm=0.3
wde=0.7
w1=-2./3.
w2=-1./2.
a0=1e-10
a1=1.0
n=1000

def f(a,t,w):
    return 1/(H*(wm/a+wr/a**2+wde/a**(3*w+1))**.5)

h=(a1-a0)/(n-1)
a=linspace(a0,a1,n)
t1=zeros(n)
t2=zeros(n)
for i in range(n-1):
    t1[i+1]=t1[i]+h*f(a[i]+.5*h,t1[i]+.5*h*f(a[i],t1[i],w1),w1)
    t2[i+1]=t2[i]+h*f(a[i]+.5*h,t2[i]+.5*h*f(a[i],t2[i],w2),w2)

plot(t1,a,label="w=-2/3")
plot(t2,a,label="w=-1/2")
xlabel("t [seconds]")
ylabel("a(t)")
legend()
show()
