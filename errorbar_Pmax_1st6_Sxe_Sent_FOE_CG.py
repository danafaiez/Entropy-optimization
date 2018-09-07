#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy.special import factorial
#Probability of localizing 3 particles in the first 6 sites is maximized using unitary.c, 100 times, each time starting with a different random set of complex gaussian coefficients in thermal_psi
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, seed1, complex gaussian coeff, x_init=0,x_fin=5, sob=6, #bath sites=12. the corresponsing S_xe, Sent, and FOE are then found at each iteration. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

P = [None]*13
Sxe = [None]*5
FOE = [None]*5
Sent = [None]*5

P[1] = 'Pmax_L10_Ne3_CG_loc1st6_B0.01.d'
P[2] = 'Pmax_L12_Ne3_CG_loc1st6_B0.01.d'
P[3] = 'Pmax_L14_Ne3_CG_loc1st6_B0.01.d'
P[4] = 'Pmax_L16_Ne3_CG_loc1st6_B0.01.d'
P[5] = 'Pmax_L18_Ne3_CG_loc1st6_B0.01.d'
P[6] = 'Pmax_L20_Ne3_CG_loc1st6_B0.01.d'
P[7] = 'Pmax_L22_Ne3_CG_loc1st6_B0.01.d'
P[8] = 'Pmax_L24_Ne3_CG_loc1st6_B0.01.d'
P[9] = 'Pmax_L26_Ne3_CG_loc1st6_B0.01.d'
P[10] = 'Pmax_L28_Ne3_CG_loc1st6_B0.01.d'
#P[11] = 'Pmax_L30_Ne3_CG_loc1st6_B0.01.d' #run this again and use it instead of next line
P[11] = 'SMALL_Pmax_L30_Ne3_CG_loc1st6_B0.01.d'

#Sxe[1] = 'Sxe_L12_Ne3_CG_loc1st6_B0.01.d'
Sxe[1] = 'Sxe_L18_Ne3_CG_loc1st6_B0.01.d'
Sxe[2] = 'Sxe_L24_Ne3_CG_loc1st6_B0.01.d'
Sxe[3] = 'Sxe_L30_Ne3_CG_loc1st6_B0.01.d'#chnage this-this is only ran for a small number of times


#Sent[1] = 'Sent_L12_Ne3_CG_loc1st6_B0.01.d'#I dont think this is done correctly since Sent find entropy btw the 1st #bath_sites and the rest where #bath_sites>6. and this wouldnt work in the case where particlesare localized in the 1st 6 sites.
Sent[1] = 'Sent_L18_Ne3_CG_loc1st6_B0.01.d'
Sent[2] = 'Sent_L24_Ne3_CG_loc1st6_B0.01.d'
Sent[3] = 'Sent_L30_Ne3_CG_loc1st6_B0.01.d'

#FOE[1] = 'FOE_L12_Ne3_CG_loc1st6_B0.01.d'
FOE[1] = 'FOE_L18_Ne3_CG_loc1st6_B0.01.d'
FOE[2] = 'FOE_L24_Ne3_CG_loc1st6_B0.01.d'
FOE[3] = 'FOE_L30_Ne3_CG_loc1st6_B0.01.d'#chnage this-this is only ran for a small number of times

#if L=12 is included:
#Sxe_ave = [4.98781405,6.2788,7.1850,7.8851]
#FOE_ave = [4.9868679,6.2878021,7.1913005,7.88574335] 
#Sent_ave= [2.47113735,2.43059225,2.15467225,1.8890501]
#if L=12 is not included
Sxe_ave = [6.2788,7.1850,7.8851]
FOE_ave = [6.2878021,7.1913005,7.88574335] 
Sent_ave= [2.43059225,2.15467225,1.8890501]

#Save_L30=7.8851
#Save_L24=7.1850
#Save_L18=6.2788
#Save_L12=4.9878

#Sent_L30=1.8890501
#Sent_L24=2.1546722
#Sent_L18=2.43059225
#Sent_L12=2.47113735

#FOE_ave_L30=7.8857433
#FOE_ave_L24=7.1913005
#FOE_ave_L18=6.2878021
#FOE_ave_L12=4.9868679

L=np.array([18,24,30])
l=6
Ne=3
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2=np.divide(N,np.power(M,2))

Ls=np.array([10,12,14,16,18,20,22,24,26,28,30])
Ns=factorial(Ls)/(factorial(Ne)*factorial(Ls-Ne))
Ms=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2s=np.divide(Ns,np.power(Ms,2))

for i in range(1,12):     
       a = loadtxt(P[i])
       plt.errorbar(NM2s[i-1],mean(a), yerr=np.array([[mean(a)-min(a) ,max(a)-mean(a)]]).T,fmt='g.',elinewidth=0.9,
                ms=4,capsize=1,label='Pmax'if i == 1 else "")


for i in range(1,4):
       a = loadtxt(Sxe[i])
       R_mean= np.divide(mean(a),Sxe_ave[i-1])
       R_min= np.divide(min(a),Sxe_ave[i-1])
       R_max= np.divide(max(a),Sxe_ave[i-1])
      
       plt.errorbar(NM2[i-1],R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='rs',elinewidth=0.9,
                ms=3,capsize=2,label='R=Sxe(Pmax)/S(xE,ave)'if i == 1 else "")

for i in range(1,4):
       a = loadtxt(FOE[i])
       R_mean= np.divide(mean(a),FOE_ave[i-1])
       R_min= np.divide(min(a),FOE_ave[i-1])
       R_max= np.divide(max(a),FOE_ave[i-1])

       plt.errorbar(NM2[i-1],R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='ys',elinewidth=0.9,
               ms=3,capsize=2,label='R=FOE(Pmax)/FOE(ave)'if i == 1 else "")



for i in range(1,4):
       a = loadtxt(Sent[i])
       R_mean= np.divide(mean(a),Sent_ave[i-1])
       R_min= np.divide(min(a),Sent_ave[i-1])
       R_max= np.divide(max(a),Sent_ave[i-1])

       plt.errorbar(NM2[i-1],R_mean,yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='ms',elinewidth=0.9,
               ms=3,capsize=2,label='R=Sent(Pmax)/Sent(ave)\nL-bath sites=6'if i == 1 else "")



#Ticks
ax.set_axisbelow(True)
ax.minorticks_on()
ax.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax.tick_params(which='both',top='on',left='on',right='on', bottom='on') 


# Adding plotting parameters
plt.legend()
#plt.title('P maxed in middle 6 sites_cmplx Gaussian coeff__N=3_B=0.01')
ax.set_xlabel('$\\frac{N}{M^{2}}$', fontsize=12)
ax.set_yticks(np.arange(0.5,1.05,0.1))
ax.set_xticks(np.arange(0, max(NM2), 3))
show()
