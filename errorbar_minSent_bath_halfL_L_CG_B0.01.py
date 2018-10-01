import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import factorial
from scipy.optimize import curve_fit
import sympy as sym
#Minimizing S_ent with #bath_sites=L/2 for different system sized, for 20 times 
#for 3 times, each time starting from a different random pure state and going through the minimization process 100 times, each time 
#starting with a differnt set of phases to do the minimization process in min.c, in unitary file. 
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, seed1, Ne=2, complex gaussian coeff. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

Sent = [None]*5
min_Sent = [None]*5
Pmax = [None]*5

min_Sent[1] = 'minSent_bath_halfL_L8_CG_B0.01.d'
min_Sent[2] = 'minSent_bath_halfL_L12_CG_B0.01.d'
min_Sent[3] = 'minSent_bath_halfL_L16_CG_B0.01.d'
min_Sent[4] = 'minSent_bath_halfL_L20_CG_B0.01.d'

Sent_ave=[1.482141,1.71113,1.850381,1.981776]
"""
Sent_ave[8]   = 1.482141
Sent_ave[L12] = 1.71113
Sent_ave[L16] = 1.850381 
Sent_ave[L20] = 1.981776
"""


Pmax[1] = 'Pmax_L8_loc4.CG.d'
Pmax[2] = 'Pmax_L12_loc6.CG.d'
Pmax[3] = 'Pmax_L16_loc8.CG.d'
Pmax[4] = 'Pmax_L20_loc10.CG.d'

Sent[1]= 'Sent_L8_loc4_CG.d'
Sent[2]= 'Sent_L12_loc6_CG.d'
Sent[3]= 'Sent_L16_loc8_CG.d'
Sent[4]= 'Sent_L20_loc10_CG.d'

#L=[8,12,16,20]
L=np.array([8,12,16,20])
l=np.subtract(L,4)
Ne=2
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne)) #M = Hilbert dim of bath
NM2=np.divide(N,np.power(M,2))
NM=np.divide(M,N)
R=[]

for i in range(1,5):
       a = loadtxt(min_Sent[i])
       R_mean= np.divide(mean(a),Sent_ave[i-1])
       R_min= np.divide(min(a),Sent_ave[i-1])
       R_max= np.divide(max(a),Sent_ave[i-1])
      
       plt.errorbar(L[i-1],R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='bs',elinewidth=0.9,
                ms=3,capsize=2,label='R=Sent(min)/Sent(ave)'if i == 1 else "")

       R.append(R_mean)
for i in range(1,5):
       a = loadtxt(Pmax[i])
       plt.errorbar(L[i-1],mean(a), yerr=np.array([[mean(a)-min(a),max(a)-mean(a)]]).T,fmt='k.',elinewidth=0.6,
                ms=4,capsize=1,label='Pmax'if i ==1 else "")



for i in range(1,5):
       a = loadtxt(Sent[i])
       R_mean= np.divide(mean(a),Sent_ave[i-1])
       R_min= np.divide(min(a),Sent_ave[i-1])
       R_max= np.divide(max(a),Sent_ave[i-1])

       plt.errorbar(L[i-1],R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='ms',elinewidth=0.9, 
               ms=3,capsize=2,label='R=Sent(loc)/Sent(ave)'if i == 1 else "")


def func(NM, a,b):
    return np.exp(-NM*a)*b
popt, pcov = curve_fit(func,R ,NM)

#plt.plot(NM, func(NM, *popt),'r--',label="Fitted Curve")
#ax.text(0.35, 0.45,'fit:bexp(-ax)/x ,  a=%s , b=%s' % (popt[0] , popt[1]) ,fontsize=8)
#Ticks
ax.set_axisbelow(True)
ax.minorticks_on()
ax.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax.tick_params(which='both',top='on',left='on',right='on', bottom='on') 

# Adding plotting parameters
plt.legend()
ax.set_xlabel('M/N', fontsize=12)
ax.set_ylabel('R', fontsize=12)
ax.set_yticks(np.arange(0,1,0.1))
#ax.set_xticks(NM)
show()
