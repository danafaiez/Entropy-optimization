#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import factorial
#Minimizing S_xe with sob=4 for different system sized, for 100 times, each time starting with a different random set of complex gaussian coefficients in thermal_psi, and each time
#starting with a differnt set of phases to do the minimization process in min.c, in unitary file. 
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, seed1, complex gaussian coeff, sob=4. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

Sxe = [None]*5

Sxe[1] = 'minSxe_L8_CG_B0.01.d'
Sxe[2] = 'minSxe_L12_CG_B0.01.d'
Sxe[3] = 'minSxe_L16_CG_B0.01.d'
Sxe[4] = 'minSxe_L20_CG_B0.01.d'

Sxe_ave = [2.9369147,3.8018902,4.3751015,4.8310735]

"""
Save_dom[8]   = 2.9369147
Save_dom[L12] = 3.8018902
Save_dom[L16] = 4.3751015
Save_dom[L20] = 4.8310735 
Save_dom[L24] = 5.2016302 #check this one
"""

L=[8,12,16,20]


for i in range(1,5):
       a = loadtxt(Sxe[i])
       R_mean= np.divide(mean(a),Sxe_ave[i-1])
       R_min= np.divide(min(a),Sxe_ave[i-1])
       R_max= np.divide(max(a),Sxe_ave[i-1])
      
       plt.errorbar(L[i-1],R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='rs',elinewidth=0.9,
                ms=3,capsize=2,label='R=Sxe(Pmax)/S(xE,ave)'if i == 1 else "")


#Ticks
ax.set_axisbelow(True)
ax.minorticks_on()
ax.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax.tick_params(which='both',top='on',left='on',right='on', bottom='on') 


# Adding plotting parameters
plt.legend()
#plt.title('Sxe is minimized with sob=4, cmplx Gaussian coeff, 2 particles, B=0.01')
ax.set_xlabel('L', fontsize=12)
ax.set_yticks(np.arange(0.5,1.05,0.1))
ax.set_xticks(L)
show()
