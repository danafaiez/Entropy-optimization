#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import factorial
#with NO Hopping allowed between the subsystem and bath, minimizing S_ent with L-subsys=4 for different system sized, for 6 times, each time starting from a different random 
#pure state and going through the minimization process 50 times, each time 
#starting with a differnt set of phases to do the minimization process in min.c, in unitary file. 
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, seed1, Ne=2, complex gaussian coeff.

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

Sent_max= [None]*5
Sent_min= [None]*5

#6 max values of Sent, each from a different initial state, ran for t=10000, delta_t=5
Sent_max[0]  = 'NH_maxSent_L8_CG_B0.01.d'#done
Sent_max[1] = 'NH_maxSent_L12_CG_B0.01.d'#done
Sent_max[2] = 'NH_maxSent_L16_CG_B0.01.d'#done
Sent_max[3] = 'NH_maxSent_L20_CG_B0.01.d'#done
Sent_max[4] = 'NH_maxSent_L24_CG_B0.01.d'#done


#6 min values for S_ent, each from a different initial state and each minimization is done for 50 times and the min of the 50 minima is found
Sent_min[0] = 'NH_new_minSent_L8_CG_B0.01.d'#done
Sent_min[1] = 'NH_new_minSent_L12_CG_B0.01.d'#done
Sent_min[2] = 'NH_new_minSent_L16_CG_B0.01.d'#done
Sent_min[3] = 'NH_new_minSent_L20_CG_B0.01.d'#done
Sent_min[4] = 'NH_new_minSent_L24_CG_B0.01.d'#done :starting from 6 different thermal states,but only doing the minimization for 20 times as supposed to 50.



L=np.array([8,12,16,20,24])
#plotting the min only
for i in range(0,5):
        min = loadtxt(Sent_min[i])
        Rent_mid=np.mean(min)
        ent_std=[]
        ent_std.append(np.std(min))
        plt.errorbar(L[i],Rent_mid,yerr=np.array([ent_std,ent_std]).T,fmt='*',color='k',elinewidth=0.6,ms=4,capsize=2,label=r'$\mathrm{min(S_{ent})};\Delta=4$' if i == 1 else "")
#plotting the max only
for i in range(0,5):
        max = loadtxt(Sent_max[i])
        Rent_mid_max=np.mean(max)
        ent_stdmax=[]
        ent_stdmax.append(np.std(max))
        plt.errorbar(L[i],Rent_mid_max,yerr=np.array([ent_stdmax,ent_stdmax]).T,fmt='*',color='orange',elinewidth=0.6,ms=4,capsize=2,label=r'$\mathrm{max(S_{ent})};\Delta=4$' if i == 1 else "")
#plotting the ratio
for i in range(0,5):
        max = loadtxt(Sent_max[i])
        min = loadtxt(Sent_min[i])
        Rent=[]
        ent_std=[]
        Rent.append(np.divide(min,max)) 
        Rent_mid=np.mean(Rent)
        ent_std.append(np.std(Rent))     
        plt.errorbar(L[i],Rent_mid,yerr=np.array([ent_std,ent_std]).T,fmt='d',color='lightcoral',elinewidth=0.6,ms=2,capsize=2,label=r'$\mathrm{R(S_{ent})};\Delta=4$'if i == 1 else "")
plt.legend()
ax.set_xlabel('L', fontsize=15)
ax.set_ylabel('R', fontsize=15)
ax.tick_params(direction='out', colors='k', labelsize=14)
ax.set_xticks(np.arange(8,25,1))
plt.savefig("NH_min_Sent_Sxe.eps",format='eps',bbox_inches='tight')
show()
