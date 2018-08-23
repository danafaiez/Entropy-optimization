#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
# Importing colormap module
import matplotlib.cm as cm
#Probability of localizing 3 particles in the middle 6 sites is maximized using unitary.c, 100 times, each time starting with a different random set of complex gaussian coefficients
#with parameters: -t 1.90 -tp 1.90 -U 0.50 -Up 0.50, B=0.01, complex gaussian coeff, sob=bath sites=6. the corresponsing S_xe is then found at each iteration. 


fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)

Y = [None]*4
S = [None]*4
Y[1] = np.loadtxt('Pmax_L30_Ne3_bath6_complG_locmid_B0.01.d')
Y[2] = np.loadtxt('Pmax_L24_Ne3_bath6_complG_locmid_B0.01.d')
Y[3] = np.loadtxt('Pmax_L18_Ne3_bath6_complG_locmid_B0.01.d')

S[1] = np.loadtxt('Smax_L30_Ne3_bath6_complG_locmid_B0.01.d')
S[2] = np.loadtxt('Smax_L24_Ne3_bath6_complG_locmid_B0.01.d')
S[3] = np.loadtxt('Smax_L18_Ne3_bath6_complG_locmid_B0.01.d')

Save= [7.8851, 7.1850, 6.2788]

#Save_L30=7.8851
#Save_L24=7.1850
#Save_L18=6.2788



# Appending RGBA values to a new array
#for i in range(20):
#    color_array.append(cm.get_cmap('tab20')(i))


for i in range(1,4):     
       a = loadtxt(Y[i])
       plt.errorbar(i,mean(a), yerr=np.array([[mean(a)-min(a) ,max(a)-mean(a)]]).T,fmt='g.',elinewidth=0.2,
                ms=4,capsize=1,label='Pmax'if i == 1 else "")


for i in range(1,4):
       a = loadtxt(S[i])
       R_mean= np.divide(Save[i-1], mean(a))
       R_min= np.divide(Save[i-1], min(a))
       R_max= np.divide(Save[i-1], max(a))
      
       plt.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min ,R_max-R_mean]]).T,fmt='rs',elinewidth=0.2,
                ms=3,capsize=2,label='R=S(xE,ave)/S(xE,min)'if i == 1 else "")


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
x=[1,2,3]
x_label=[18,24,30]
ax.set_xlabel('L', fontsize=12)
y_label=np.arange(0.5,2.5,0.5)
ax.set_yticks(y_label)
ax.set_xticks(x)
ax.set_xticklabels(x_label)
show()
