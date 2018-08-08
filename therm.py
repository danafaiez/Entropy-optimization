#!/Library/Frameworks/Python.framework/Versions/Current/bin/python
from scipy import *
from pylab import *
import scipy.optimize
import sys

def test_E(n):
   e = zeros((n))
   for i in range(n):
      e[i] = bin(i).count("1")
   return sort(e)
#print(e(1))
def plot_test(s_out,n):
   plot(s_out[:,0],s_out[:,1])
   x = s_out[:,0]
   y = x/n
   plot(x,-n*(y*log(y)+(1-y)*log(1-y)))
   draw()
   show()

def Z(beta):
   return sum(exp(-beta*e))
def F(beta):
   return -log(Z(beta))/beta

def E(beta):
   return sum(exp(-beta*e)*e)/Z(beta)
#F=E-TS, 
def S(beta):
   return beta*(E(beta)-F(beta))
def f(beta):
   return energy -E(beta)

num_skip=2

test_m = 0
if test_m:
   e = test_E(1<<test_m)
else:
   e = loadtxt(sys.argv[1])
N = len(e)
m = int(N/num_skip)
s_out = zeros((m,3))

e2 = reshape(e[0:m*num_skip],(m,num_skip))
e_ave = mean(e2,axis=1)

i=0
for energy in e_ave:
   b = scipy.optimize.bisect(f,-30,30)
   s_out[i,0] = energy
   s_out[i,1] = S(b)
   s_out[i,2] = b
   i += 1
if test_m:
    plot_test(s_out,test_m)
    savetxt("S_th_test.dat",s_out)
else:
    savetxt("S_th.dat",s_out)


