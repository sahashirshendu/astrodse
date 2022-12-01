from pylab import *

H=70/3.086e19
wr=5e-5
wm=0.3
wk1=-0.3
wk2=-0.5
wk3=-0.7
a0=1e-10
a1=1.0
n=1000

def f(a,t,wk):
    return 1/(H**2*(wm/a+wr/a**2-wk))**.5

h=(a1-a0)/(n-1)
a=linspace(a0,a1,n)
t1=zeros(n)
t2=zeros(n)
t3=zeros(n)
for i in range(n-1):
    t1[i+1]=t1[i]+h*f(a[i]+.5*h,t1[i]+.5*h*f(a[i],t1[i],wk1),wk1)
    t2[i+1]=t2[i]+h*f(a[i]+.5*h,t2[i]+.5*h*f(a[i],t2[i],wk2),wk2)
    t3[i+1]=t3[i]+h*f(a[i]+.5*h,t3[i]+.5*h*f(a[i],t3[i],wk3),wk3)

plot(t1,a,'-',label='$\\Omega_k = -0.3$')
plot(t2,a,'-',label='$\\Omega_k = -0.5$')
plot(t3,a,'-',label='$\\Omega_k = -0.7$')
xlabel('t [seconds]')
ylabel('a(t)')
legend()
show()
