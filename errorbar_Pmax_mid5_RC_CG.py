#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
from scipy.special import factorial
import matplotlib.cm as cm
#Probability of localizing 3 particles in the middle 5 sites is maximized using unitary.c, 100 times, each time starting with a different random set of real/complex gaussian coefficients
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, real gaussian coeff, sob=bath sites=6. the corresponsing S_xe is then found at each iteration. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

PR = [None]*10
PC = [None]*10
PR[9] = np.loadtxt('Pmax_L30_Ne3_RG_locmid5_B0.01.d')
PR[8] = np.loadtxt('Pmax_L27_Ne3_RG_locmid5_B0.01.d')
PR[7] = np.loadtxt('Pmax_L25_Ne3_RG_locmid5_B0.01.d')
PR[6] = np.loadtxt('Pmax_L22_Ne3_RG_locmid5_B0.01.d')
PR[5] = np.loadtxt('Pmax_L20_Ne3_RG_locmid5_B0.01.d')
PR[4] = np.loadtxt('Pmax_L17_Ne3_RG_locmid5_B0.01.d')
PR[3] = np.loadtxt('Pmax_L15_Ne3_RG_locmid5_B0.01.d')
PR[2] = np.loadtxt('Pmax_L12_Ne3_RG_locmid5_B0.01.d')
PR[1] = np.loadtxt('Pmax_L10_Ne3_RG_locmid5_B0.01.d')

PC[9] = np.loadtxt('Pmax_L30_Ne3_CG_locmid5_B0.01.d')
PC[8] = np.loadtxt('Pmax_L27_Ne3_CG_locmid5_B0.01.d')
PC[7] = np.loadtxt('Pmax_L25_Ne3_CG_locmid5_B0.01.d')
PC[6] = np.loadtxt('Pmax_L22_Ne3_CG_locmid5_B0.01.d')
PC[5] = np.loadtxt('Pmax_L20_Ne3_CG_locmid5_B0.01.d')
PC[4] = np.loadtxt('Pmax_L17_Ne3_CG_locmid5_B0.01.d')
PC[3] = np.loadtxt('Pmax_L15_Ne3_CG_locmid5_B0.01.d')
PC[2] = np.loadtxt('Pmax_L12_Ne3_CG_locmid5_B0.01.d')
PC[1] = np.loadtxt('Pmax_L10_Ne3_CG_locmid5_B0.01.d')


L=np.array([10,12,15,17,20,22,25,27,30])
l=5
Ne=3
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2=np.divide(N,np.power(M,2))

for i in range(1,10):     
       a = loadtxt(PR[i])
       plt.errorbar(NM2[i-1],mean(a), yerr=np.array([[mean(a)-min(a) ,max(a)-mean(a)]]).T,color='crimson',fmt='.',elinewidth=0.6,
                ms=4,capsize=1,label='Real Gaussian coeff'if i == 1 else "")

for i in range(1,10):
       a = loadtxt(PC[i])
       plt.errorbar(NM2[i-1],mean(a), yerr=np.array([[mean(a)-min(a) ,max(a)-mean(a)]]).T,fmt='k.',elinewidth=0.6,
                ms=4,capsize=1,label='Complex Gaussian coeff'if i == 1 else "")

#Ticks
ax.set_axisbelow(True)
ax.minorticks_on()
ax.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')
ax.tick_params(which='both',top='off',left='on',right='on', bottom='on') 

#ax.legend(loc='best', frameon=True)

# Adding plotting parameters
plt.legend()
#plt.title('P maxed in middle 5 sites with real and complex Gaussian coeff;N=3;B=0.01',fontsize=10)
ax.set_xlabel('$\\frac{N}{M^{2}}$', fontsize=12)
ax.set_ylabel('$P_{max}$', fontsize=12)
ax.set_yticks(np.arange(0.4,1.05,0.1))
ax.set_xticks(np.arange(0, np.add(max(NM2),3), 3))
plt.savefig("error_Pmax_mid5.png",bbox_inches='tight')
show()
