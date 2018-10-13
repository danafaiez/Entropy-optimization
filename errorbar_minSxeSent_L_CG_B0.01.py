#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import factorial
#Minimizing S_xe with sob=4,L/2 for different system sized, for 3 times, each time starting from a different random pure state and going through the minimization process 100 times, each time 
#Doing the same for S_ent with L-bath=4 or L/2.
#starting with a differnt set of phases to do the minimization process in min.c, in unitary file. 
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, seed1, Ne=2, complex gaussian coeff, sob=4. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

Sxe = [None]*5
Sent = [None]*5
Sent_half= [None]*5
Sxe[1] = 'minSxe_L8_CG_B0.01.d'
Sxe[2] = 'minSxe_L12_CG_B0.01.d'
Sxe[3] = 'minSxe_L16_CG_B0.01.d'
Sxe[4] = 'minSxe_L20_CG_B0.01.d'

Sxe_ave = [2.9369147,3.8018902,4.3751015,4.8310735]

"""
SxE_ave[8]   = 2.9369147
SxE_ave[L12] = 3.8018902
SxE_ave[L16] = 4.3751015
SxE_ave[L20] = 4.8310735 
SxE_ave[L24] = 5.2016302 #check this one
"""

Sent[1] = 'minSent_L8_bath4_CG_B0.01.d'
Sent[2] = 'minSent_L12_bath8_CG_B0.01.d'
Sent[3] = 'minSent_L16_bath12_CG_B0.01.d'
Sent[4] = 'minSent_L20_bath16_CG_B0.01.d'

Sent_ave = [1.4821405,1.4680688,1.32090275,1.19370775]
"""
S_th_loc4=1.002194
#corresponding Beta:0.5 (compared to B=0.01)
S_th_tot=[1.251,1.618,2.368,2.799]
#corresponding Beta=[0.54,0.568,0.506,0.501]
R_S_th=np.divide(S_th_loc4,S_th_tot)

Sent_ave[8]   = 1.4821405
Sent_ave[L12] = 1.4680688
Sent_ave[L16] = 1.32090275
Sent_ave[L20] = 1.19370775
Sent_ave[L24] =?
"""
#bath=sob=L/2
#only done for 20 times. redo this for 300 times.
Sent_half[1] = 'minSent_bath_halfL_L8_CG_B0.01.d'
Sent_half[2] = 'minSent_bath_halfL_L12_CG_B0.01.d'
Sent_half[3] = 'minSent_bath_halfL_L16_CG_B0.01.d'
Sent_half[4] = 'minSent_bath_halfL_L20_CG_B0.01.d'

Sent_half_ave=[1.482141,1.71113,1.850381,1.981776]
"""
Sent_ave[8]   = 1.482141
Sent_ave[L12] = 1.71113
Sent_ave[L16] = 1.850381 
Sent_ave[L20] = 1.981776
"""


L=np.array([8,12,16,20])
l=4
Ne=2
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2=np.divide(N,np.power(M,2))


for i in range(1,5):
       a = loadtxt(Sxe[i])
       R_mean= np.divide(mean(a),Sxe_ave[i-1])
       R_min= np.divide(min(a),Sxe_ave[i-1])
       R_max= np.divide(max(a),Sxe_ave[i-1])
      
       plt.errorbar(L[i-1],R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='g*',elinewidth=0.6,
               ms=3,capsize=2,label=r'$R(S_{xE})$,$\Delta{x}=4$'if i == 1 else "")


for i in range(1,5):
       a = loadtxt(Sent[i])
       R_mean= np.divide(mean(a),Sent_ave[i-1])
       R_min= np.divide(min(a),Sent_ave[i-1])
       R_max= np.divide(max(a),Sent_ave[i-1])

       plt.errorbar(L[i-1],R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='m*',elinewidth=0.6,
                ms=3,capsize=2,label=r'$R(S_{ent})$,Bath$=L-4$'if i == 1 else "")
"""
for i in range(1,5):
       a = loadtxt(Sent_half[i])
       R_mean= np.divide(mean(a),Sent_half_ave[i-1])
       R_min= np.divide(min(a),Sent_half_ave[i-1])
       R_max= np.divide(max(a),Sent_half_ave[i-1])

       plt.errorbar(N[i-1],R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='b*',elinewidth=0.6,
                ms=3,capsize=2,label=r'$R(S_{ent})$,Bath$=L/2$'if i == 1 else "")
"""
#plt.scatter(N,R_S_th)
# Adding plotting parameters
plt.legend()
#plt.title('Sxe is minimized with sob=4, cmplx Gaussian coeff, 2 particles, B=0.01')
ax.set_xlabel('L', fontsize=12)
ax.set_ylabel('R', fontsize=12)
#ax.set_yticks(np.arange(0.5,1.05,0.1))
ax.set_yticks(np.arange(0,1.1,0.1))
ax.set_xticks(L)
plt.savefig("min_Sent_Sxe",bbox_inches='tight')
show()
