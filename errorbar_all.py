#!usr/local.bin/python
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
# Importing colormap module
import matplotlib.cm as cm
fig, (ax1, ax2, ax3, ax4) = plt.subplots(nrows=4, ncols=1, sharex=False)

subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=1.2)

Y = [None]*4
S = [None]*4
Y[3] = 'Pmax_L30_Ne3_bath5_complG_locside_B0.01.d'
Y[2] = 'Pmax_L25_Ne3_bath5_complG_locside_B0.01.d'
Y[1] = 'Pmax_L15_Ne3_bath5_complG_locside_B0.01.d'
S[3] = 'Smin_L30_Ne3_bath5_complG_locside_B0.01.d'
S[2] = 'Smin_L25_Ne3_bath5_complG_locside_B0.01.d'
S[1] = 'Smin_L15_Ne3_bath5_complG_locside_B0.01.d'
Save= [5.697788,7.316264,7.882103]


YY = [None]*4
SS = [None]*4
YY[3] = 'Pmax_L30_Ne3_bath6_complG_locmid_B0.01.d'
YY[2] = 'Pmax_L24_Ne3_bath6_complG_locmid_B0.01.d'
YY[1] = 'Pmax_L18_Ne3_bath6_complG_locmid_B0.01.d'
SS[3] = 'Smin_L30_Ne3_bath6_complG_locmid_B0.01.d'
SS[2] = 'Smin_L24_Ne3_bath6_complG_locmid_B0.01.d'
SS[1] = 'Smin_L18_Ne3_bath6_complG_locmid_B0.01.d'
SSave= [6.2788,7.1850,7.8851]

Ym5 = [None]*4
Sm5 = [None]*4
Ym5[3] = 'Pmax_L30_Ne3_bath5_complG_locmid_B0.01.d'
Ym5[2] = 'Pmax_L25_Ne3_bath5_complG_locmid_B0.01.d'
Ym5[1] = 'Pmax_L15_Ne3_bath5_complG_locmid_B0.01.d'
Sm5[3] = 'Smin_L30_Ne3_bath5_complG_locmid_B0.01.d'
Sm5[2] = 'Smin_L25_Ne3_bath5_complG_locmid_B0.01.d'
Sm5[1] = 'Smin_L15_Ne3_bath5_complG_locmid_B0.01.d'
Sm5ave= [5.697788,7.316264,7.882103]

Y5mR = [None]*4
S5mR = [None]*4
Y5mR[3] = np.loadtxt('Pmax_L30_Ne3_bath5_realG_locmid_B0.01.d')
Y5mR[2] = np.loadtxt('Pmax_L25_Ne3_bath5_realG_locmid_B0.01.d')
Y5mR[1] = np.loadtxt('Pmax_L15_Ne3_bath5_realG_locmid_B0.01.d')
S5mR[3] = np.loadtxt('Smin_L30_Ne3_bath5_realG_locmid_B0.01.d')
S5mR[2] = np.loadtxt('Smin_L25_Ne3_bath5_realG_locmid_B0.01.d')
S5mR[1] = np.loadtxt('Smin_L15_Ne3_bath5_realG_locmid_B0.01.d')
S5mRave= [5.6554212,7.2866015,7.8575661]

for i in range(1,4):
       a = loadtxt(Y[i])
       ax1.errorbar(i,mean(a), yerr=np.array([[mean(a)-min(a),max(a)-mean(a)]]).T,fmt='g.',elinewidth=0.2,ms=4,capsize=1,label='Pmax'if i == 1 else "")


for i in range(1,4):
       a = loadtxt(S[i])
       R_mean= np.divide(Save[i-1], mean(a))
       R_min= np.divide(Save[i-1], min(a))
       R_max= np.divide(Save[i-1], max(a))
       ax1.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='rs',elinewidth=0.9,ms=3,capsize=2,label='R=S(xE,ave)/S(xE,min)'if i == 1 else "")

ax1.set_axisbelow(True)
ax1.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax1.set_title('P maxed in first 5 sites_cmplx Gaussian coeff__N=3_B=0.01',fontsize=6)
x=[1,2,3]
x_label=[15,25,30]
ax1.set_xlabel('L', fontsize=12)
y_label=np.arange(0.5,2.5,0.5)
ax1.set_yticks(y_label)
ax1.set_xticks(x)
ax1.set_xticklabels(x_label)

for i in range(1,4):
       a = loadtxt(YY[i])
       ax2.errorbar(i,mean(a), yerr=np.array([[mean(a)-min(a) ,max(a)-mean(a)]]).T,fmt='g.',elinewidth=0.2,ms=4,capsize=1,label='Pmax'if i == 1 else "")

for i in range(1,4):
       a = loadtxt(SS[i])
       R_mean= np.divide(SSave[i-1], mean(a))
       R_min= np.divide(SSave[i-1], min(a))
       R_max= np.divide(SSave[i-1], max(a))
       ax2.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='rs',elinewidth=0.9,ms=3,capsize=2,label='R=S(xE,ave)/S(xE,min)'if i == 1 else "")

ax2.set_axisbelow(True)
ax2.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax2.set_title('P maxed in middle 6 sites_cmplx Gaussian coeff__N=3_B=0.01',fontsize=6)
x=[1,2,3]
x_label=[18,24,30]
ax2.set_xlabel('L', fontsize=12)
y_label=np.arange(0.5,2.5,0.5)
ax2.set_yticks(y_label)
ax2.set_xticks(x)
ax2.set_xticklabels(x_label)

for i in range(1,4):
       a = loadtxt(Ym5[i])
       ax3.errorbar(i,mean(a), yerr=np.array([[mean(a)-min(a),max(a)-mean(a)]]).T,fmt='g.',elinewidth=0.2,ms=4,capsize=1,label='Pmax'if i == 1 else "")


for i in range(1,4):
       a = loadtxt(Sm5[i])
       R_mean= np.divide(Sm5ave[i-1], mean(a))
       R_min= np.divide(Sm5ave[i-1], min(a))
       R_max= np.divide(Sm5ave[i-1], max(a))
       ax3.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='rs',elinewidth=0.9,ms=3,capsize=2,label='R=S(xE,ave)/S(xE,min)'if i == 1 else "")

ax3.set_axisbelow(True)
ax3.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax3.set_title('P maxed in middle 5 sites_cmplx Gaussian coeff__N=3_B=0.01',fontsize=6)
x=[1,2,3]
x_label=[15,25,30]
ax3.set_xlabel('L', fontsize=12)
y_label=np.arange(0.5,2.5,0.5)
ax3.set_yticks(y_label)
ax3.set_xticks(x)
ax3.set_xticklabels(x_label)


for i in range(1,4):
       a = loadtxt(Y5mR[i])
       ax4.errorbar(i,mean(a), yerr=np.array([[mean(a)-min(a),max(a)-mean(a)]]).T,fmt='g.',elinewidth=0.2,ms=4,capsize=1,label='Pmax'if i == 1 else "")


for i in range(1,4):
       a = loadtxt(S5mR[i])
       R_mean= np.divide(S5mRave[i-1], mean(a))
       R_min= np.divide(S5mRave[i-1], min(a))
       R_max= np.divide(S5mRave[i-1], max(a))
       ax4.errorbar(i,R_mean, yerr=np.array([[R_mean-R_min,R_max-R_mean]]).T,fmt='rs',elinewidth=0.2, ms=3,capsize=2,label='R=S(xE,ave)/S(xE,min)'if i == 1 else "")
ax4.set_axisbelow(True)
ax4.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
ax4.set_title('P maxed in middle 5 sites_real Gaussian coeff__N=3_B=0.01',fontsize=6)
x=[1,2,3]
x_label=[15,25,30]
ax4.set_xlabel('L', fontsize=12)
y_label=np.arange(0.5,2.5,0.5)
ax4.set_yticks(y_label)
ax4.set_xticks(x)
ax4.set_xticklabels(x_label)


plt.show()
