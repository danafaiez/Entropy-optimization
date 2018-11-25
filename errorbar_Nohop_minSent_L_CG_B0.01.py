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

Sxe_max= [None]*5
Sxe_min= [None]*5
"""
#6 max values of Sxe, each from a different initial state, ran for t=200, delta_t=5
Sxe_max[0]  = 'NH_maxSxe_L8_CG_B0.01.d'
Sxe_max[1] = 'NH_maxSxe_L12_CG_B0.01.d'
Sxe_max[2] = 'NH_maxSxe_L16_CG_B0.01.d'
Sxe_max[3] = 'NH_maxSxe_L20_CG_B0.01.d'

#6 min values for S_xe, each from a different initial state and each minimization is done for 50 times and the min of the 50 minima is found
Sxe_min[0] = 'NH_new_minSxe_L8_CG_B0.01.d'
Sxe_min[1] = 'NH_new_minSxe_L12_CG_B0.01.d'
Sxe_min[2] = 'NH_new_minSxe_L16_CG_B0.01.d'
Sxe_min[3] = 'NH_new_minSxe_L20_CG_B0.01.d'
"""
Sent_max= [None]*5
Sent_min= [None]*5

#6 max values of Sent, each from a different initial state, ran for t=10000, delta_t=5
Sent_max[0]  = 'NH_maxSent_L8_CG_B0.01.d'#done
Sent_max[1] = 'NH_maxSent_L12_CG_B0.01.d'#done
Sent_max[2] = 'NH_maxSent_L16_CG_B0.01.d'#done
Sent_max[3] = 'NH_maxSent_L20_CG_B0.01.d'#done

#6 min values for S_ent, each from a different initial state and each minimization is done for 50 times and the min of the 50 minima is found
Sent_min[0] = 'NH_new_minSent_L8_CG_B0.01.d'#done
Sent_min[1] = 'NH_new_minSent_L12_CG_B0.01.d'#done
Sent_min[2] = 'NH_new_minSent_L16_CG_B0.01.d'#done
Sent_min[3] = 'NH_new_minSent_L20_CG_B0.01.d'#done

"""
###below, if the same as above expect for L=6,9,12,15,18 and deta_x=3, Ne=2
#Sxe_delta3_ave= [None]*5
Sxe_delta3_max= [None]*5
Sxe_delta3_min= [None]*5

#6 max values of Sxe with delta_x=3 sites, each from a different initial state, ran for t=200, delta_t=5
Sxe_delta3_max[0] = 'NH_maxSxe_L9_delt3_CG_B0.01.d'
Sxe_delta3_max[1] = 'NH_maxSxe_L12_delt3_CG_B0.01.d'
Sxe_delta3_max[2] = 'NH_maxSxe_L15_delt3_CG_B0.01.d'
Sxe_delta3_max[3] = 'NH_maxSxe_L18_delt3_CG_B0.01.d'

#6 min values for S_xe with delta_x=3 sites, each from a different initial state and each minimization is done for 50 times and the min of the 50 minima is found
Sxe_delta3_min[0] = 'NH_new_minSxe_L9_delta3_CG_B0.01.d'
Sxe_delta3_min[1] = 'NH_new_minSxe_L12_delta3_CG_B0.01.d'
Sxe_delta3_min[2] = 'NH_new_minSxe_L15_delta3_CG_B0.01.d'
Sxe_delta3_min[3] = 'NH_new_minSxe_L18_delta3_CG_B0.01.d'

#Sent_delta3_ave= [None]*5
Sent_delta3_max= [None]*5
Sent_delta3_min= [None]*5

#6 max values of Sent for subsystem of size 3 sites, each from a different initial state, ran for t=1000, delta_t=5
Sent_delta3_max[0]  = 'NH_maxSent_L9_delta3_CG_B0.01.d'
Sent_delta3_max[1] = 'NH_maxSent_L12_delta3_CG_B0.01.d'
Sent_delta3_max[2] = 'NH_maxSent_L15_delta3_CG_B0.01.d'
Sent_delta3_max[3] = 'NH_maxSent_L18_delta3_CG_B0.01.d'

#6 min values for S_ent for subsystem of size 3 sites, each from a different initial state and each minimization is done for 50 times and the min of the 50 minima is found
Sent_delta3_min[0] = 'NH_new_minSent_L9_delta3_CG_B0.01.d'
Sent_delta3_min[1] = 'NH_new_minSent_L12_delta3_CG_B0.01.d'
Sent_delta3_min[2] = 'NH_new_minSent_L15_delta3_CG_B0.01.d'
Sent_delta3_min[3] = 'NH_new_minSent_L18_delta3_CG_B0.01.d'
"""
L=np.array([8,12,16,20])
L3=np.array([9,12,15,18])
l=4
Ne=2
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2=np.divide(N,np.power(M,2))
"""
for i in range(0,4):#change this to 4 later
        max = loadtxt(Sxe_max[i])
        min = loadtxt(Sxe_min[i])
        R_max=[]
        stdmax=[]
        R_max.append(np.divide(min,max))
        R_mid_max=np.mean(R_max)
        stdmax.append(np.std(R_max))
        plt.errorbar(L[i],R_mid_max, yerr=np.array([stdmax,stdmax]).T,fmt='b*',elinewidth=0.6,ms=4,capsize=2,label=r'$\mathrm{R(S_{xE})};\Delta=4$'if i == 1 else "")

for i in range(0,4):
        max = loadtxt(Sxe_delta3_max[i])
        min = loadtxt(Sxe_delta3_min[i])
        R_max=[]
        stdmax=[]
        R_max.append(np.divide(min,max))
        R_mid_max=np.mean(R_max)
        stdmax.append(np.std(R_max))
        plt.errorbar(L3[i],R_mid_max,yerr=np.array([stdmax,stdmax]).T,fmt='d',color='royalblue',elinewidth=0.6,ms=2,capsize=2,label=r'$\mathrm{R(S_{xE})};\Delta=3$'if i == 1 else "")
"""
for i in range(0,4):
        max = loadtxt(Sent_max[i])
        min = loadtxt(Sent_min[i])
        Rent_max=[]
        ent_stdmax=[]
        Rent_max.append(np.divide(min,max)) 
        Rent_mid_max=np.mean(Rent_max) 
        ent_stdmax.append(np.std(Rent_max))
        plt.errorbar(L[i],Rent_mid_max,yerr=np.array([ent_stdmax,ent_stdmax]).T,fmt='*',color='crimson',elinewidth=0.6,ms=4,capsize=2,label=r'$\mathrm{R(S_{ent})};\Delta=4$' if i == 1 else "")

#plotting the min only
for i in range(0,4):
        min = loadtxt(Sent_min[i])
        ent_std=[]
        Rent_mid=np.mean(min)
        ent_std.append(np.std(min))
        plt.errorbar(L[i],Rent_mid,yerr=np.array([ent_std,ent_std]).T,fmt='*',color='k',elinewidth=0.6,ms=4,capsize=2,label=r'$\mathrm{min(S_{ent})};\Delta=4$' if i == 1 else "")
#plotting the max only
for i in range(0,4):
        max = loadtxt(Sent_max[i])
        Rent_max=[]
        ent_stdmax=[]
        Rent_mid_max=np.mean(max)
        ent_stdmax.append(np.std(max))
        plt.errorbar(L[i],Rent_mid_max,yerr=np.array([ent_stdmax,ent_stdmax]).T,fmt='*',color='orange',elinewidth=0.6,ms=4,capsize=2,label=r'$\mathrm{max(S_{ent})};\Delta=4$' if i == 1 else "")


"""
for i in range(0,4):
        max = loadtxt(Sent_max[i])
        min = loadtxt(Sent_min[i])
        Rent_max=[]
        ent_stdmax=[]
        Rent_max.append(np.divide(min,max)) 
        Rent_mid_max=np.mean(Rent_max)
        ent_stdmax.append(np.std(Rent_max))     
        plt.errorbar(L3[i],Rent_mid_max,yerr=np.array([ent_stdmax,ent_stdmax]).T,fmt='d',color='lightcoral',elinewidth=0.6,ms=2,capsize=2,label=r'$\mathrm{R(S_{ent})};\Delta=3$'if i == 1 else "")
"""
plt.legend()
ax.set_xlabel('L', fontsize=15)
ax.set_ylabel('R', fontsize=15)
ax.tick_params(direction='out', colors='k', labelsize=14)
#ax.set_yticks(np.arange(0,1.1,0.2))
ax.set_xticks(np.arange(8,21,1))
plt.savefig("new_min_Sent_Sxe.eps",format='eps',bbox_inches='tight')
show()
