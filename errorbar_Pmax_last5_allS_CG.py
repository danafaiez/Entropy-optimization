#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.special import factorial
import matplotlib.cm as cm
#Probability of localizing 3 particles in the last 5 sites is maximized using unitary.c, 2 times, each time starting with a different random set of complex gaussian coefficients
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, complex gaussian coeff, L-bath sites=5. the corresponsing S_ent, FOE, and SxE is then found at each iteration. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

P = [None]*6
Sent = [None]*6
FOE = [None]*6
SxE = [None]*6
P[5] = 'Pmax_L30_Ne3_complG_loclast5_B0.01.d'
P[4] = 'Pmax_L25_Ne3_complG_loclast5_B0.01.d'
P[3] = 'Pmax_L20_Ne3_complG_loclast5_B0.01.d'
P[2] = 'Pmax_L15_Ne3_complG_loclast5_B0.01.d'
P[1] = 'Pmax_L10_Ne3_complG_loclast5_B0.01.d'

#Sent:
#redo for more data points
Sent[5] = 'Sent_L30_Ne3_complG_loclast5_B0.01.d'#this was killed after the last iter due to incorret freeing. 
Sent[4] = 'Sent_L25_Ne3_complG_loclast5_B0.01.d'
Sent[3] = 'Sent_L20_Ne3_complG_loclast5_B0.01.d'
Sent[2] = 'Sent_L15_Ne3_complG_loclast5_B0.01.d'
Sent[1] = 'Sent_L10_Ne3_complG_loclast5_B0.01.d'
Sent_ave= [2.2795262,2.283893,2.0295868,1.788395,1.597542]
#Sent_ave_L30= 1.597542 
#Sent_ave_L25= 1.788395
#Sent_ave_L20=2.0295868 
#Sent_ave_L15= 2.283893 
#Sent_ave_L10=2.2795262 

#SxE:
SxE[5] = 'SxE_L30_Ne3_complG_sob5_B0.01.d'#this was killed after the last iter due to incorret freeing. 
SxE[4] = 'SxE_L25_Ne3_complG_sob5_B0.01.d'
SxE[3] = 'SxE_L20_Ne3_complG_sob5_B0.01.d'
SxE[2] = 'SxE_L15_Ne3_complG_sob5_B0.01.d'
SxE[1] = 'SxE_L10_Ne3_complG_sob5_B0.01.d'
SxE_ave= [4.3828965,5.6948083,6.6200108,7.3174784,7.8856833]
#SxE_ave_L30=7.8856833 
#Sxe_ave_L25=7.3174784
#Sxe_ave_L20=6.6200108 
#SxE_ave_L15=5.6948083 
#SxE_ave_L10=4.3828965 

#FOE:
FOE[5] = 'FOE_L30_Ne3_complG_sob5_B0.01.d' #done remotely-killed after 15th iter
FOE[4] = 'FOE_L25_Ne3_complG_sob5_B0.01.d'
FOE[3] = 'FOE_L20_Ne3_complG_sob5_B0.01.d'
FOE[2] = 'FOE_L15_Ne3_complG_sob5_B0.01.d'
FOE[1] = 'FOE_L10_Ne3_complG_sob5_B0.01.d'
FOE_ave= [4.386126,5.690196,6.620477,7.318672,7.888582]
#FOE_ave_L30=7.888582 
#FOE_ave_L25=7.318672 
#FOE_ave_L20=6.620477 
#FOE_ave_L15=5.690196 
#FOE_ave_L10=4.386126 


L=np.array([10,15,20,25,30])
l=5
Ne=3
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2=np.divide(N,np.power(M,2))

for i in range(1,6):     
       a = loadtxt(P[i])
       plt.errorbar(i,mean(a), yerr=np.array([[mean(a)-min(a) ,max(a)-mean(a)]]).T,fmt='k.',elinewidth=0.6,
                ms=4,capsize=1,label='Pmax'if i ==1 else "")


for i in range(1,6):
       a = loadtxt(Sent[i])
       R_mean= np.divide(mean(a),Sent_ave[i-1])
       R_min= np.divide(min(a),Sent_ave[i-1])
       R_max= np.divide(max(a),Sent_ave[i-1])
      
       plt.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='m*',elinewidth=0.6,
                ms=3,capsize=2,label='R=Sent(loc)/Sent(ave)'if i == 1 else "")


for i in range(1,6):
       a = loadtxt(SxE[i])
       R_mean= np.divide(mean(a),SxE_ave[i-1])
       R_min= np.divide(min(a),SxE_ave[i-1])
       R_max= np.divide(max(a),SxE_ave[i-1])

       plt.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='g*',elinewidth=0.6,
                ms=3,capsize=2,label='R=SxE(loc)/SxE(ave)'if i == 1 else "")

for i in range(1,6):
       a = loadtxt(FOE[i])
       R_mean= np.divide(mean(a),FOE_ave[i-1])
       R_min= np.divide(min(a),FOE_ave[i-1])
       R_max= np.divide(max(a),FOE_ave[i-1])

       plt.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='r*',elinewidth=0.6,
                ms=3,capsize=2,label='R=FOE(loc)/FOE(ave)'if i == 1 else "")


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
plt.title('P maxed in last 5 sites_cmplx Gaussian coeff__N=3_B=0.01')
x=[1,2,3,4,5]
x_label=[10,15,20,25,30]
ax.set_xlabel('L', fontsize=12)
y_label=np.arange(0.4,0.85,0.05)
ax.set_yticks(y_label)
ax.set_xticks(x)
ax.set_xticklabels(x_label)
show()
