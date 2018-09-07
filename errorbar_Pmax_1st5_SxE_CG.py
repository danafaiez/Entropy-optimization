#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
# Importing colormap module
import matplotlib.cm as cm
#Probability of localizing 3 particles in the first 5 sites is maximized using unitary.c, 100 times, each time starting with a different random set of complex gaussian coefficients
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, complex gaussian coeff, sob=bath sites=5. the corresponsing S_xe is then found at each iteration. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

Y = [None]*4
S = [None]*4

Y[3] = 'Pmax_L30_Ne3_bath5_complG_locside_B0.01.d'
Y[2] = 'Pmax_L25_Ne3_bath5_complG_locside_B0.01.d'
Y[1] = 'Pmax_L15_Ne3_bath5_complG_locside_B0.01.d'

S[3] = 'Smin_L30_Ne3_bath5_complG_locside_B0.01.d'
S[2] = 'Smin_L25_Ne3_bath5_complG_locside_B0.01.d'
S[1] = 'Smin_L15_Ne3_bath5_complG_locside_B0.01.d'

Save= [5.697788,7.316264,7.882103]

#Save_L30=7.882103
#Save_L25=7.316264
#Save_L15=5.697788


for i in range(1,4):     
       a = loadtxt(Y[i])
       plt.errorbar(i,mean(a), yerr=np.array([[mean(a)-min(a) ,max(a)-mean(a)]]).T,fmt='g.',elinewidth=0.2,
                ms=4,capsize=1,label='Pmax'if i == 1 else "")


for i in range(1,4):
       a = loadtxt(S[i])
       R_mean= np.divide(mean(a),Save[i-1])
       R_min= np.divide(min(a),Save[i-1])
       R_max= np.divide(max(a),Save[i-1])
      
       plt.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='rs',elinewidth=0.9,
                ms=3,capsize=2,label='R=SxE(min)/SxE(ave)'if i == 1 else "")


ax.set_axisbelow(True)
ax.minorticks_on()
# Customize the major grid
ax.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
# Customize the minor grid
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
# display of all ticks.
ax.tick_params(which='both',top='on',left='on',right='on', bottom='on') 

#ax.legend(loc='best', frameon=True)

# Adding plotting parameters
plt.legend()
plt.title('P maxed in first 5 sites_cmplx Gaussian coeff__N=3_B=0.01')
x=[1,2,3]
x_label=[15,25,30]
ax.set_xlabel('L', fontsize=12)
y_label=np.arange(0.5,0.7,0.05)
ax.set_yticks(y_label)
ax.set_xticks(x)
ax.set_xticklabels(x_label)
show()
