#generic, using newPsi unitary method to maximize the probability of localization in the first 5 or 6 sites .
#plot of maxP vs N/M^2,  B=0.1. With Ne=3.
#!usr/local.bin/python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
from scipy.special import factorial
import matplotlib.pyplot as plt
from pylab import *
L=np.array([6,8,10,12,14,16,18,20,22,24,26,28,30])
maxp=[0.981708,0.891795,0.701720,0.655701 ,0.607702,0.569807,0.570462,0.552288,0.532385 ,0.556465 ,0.522108,0.526340,0.506590]
maxp_U0=[0.996045 ,0.822170 ,0.765447 ,0.713670 ,0.665806 ,0.613745 ,0.617839 ,0.602222 ,0.584608 ,0.561824 ,0.562963 ,0.557483,0.528033]
maxp_real=[0.854391 ,0.739394 ,0.587737 ,0.558884 ,0.509882 ,0.474655 ,0.486877 ,0.445262 ,0.434734,0.441990 ,0.416185 ,0.415721 ,0.406627]
maxp_U0_real=[0.549035 ,0.671590 ,0.664355 ,0.592277 ,0.560391 ,0.490247 ,0.498986,0.479993,0.480977 ,0.450987,0.456395,0.444646,0.418540]
##maxP for generic system, t=1.0,tp=0.96,U=0.5,Up=0.48,complex_guass coef,B=0.1,localized in 5 sites
maxp_Ulow_5=[0.985325 ,0.816349 ,0.686237,0.671494 ,0.627503 ,0.610844 ,0.594091 ,0.568750,0.545662 ,0.550307 ,0.553379 ,0.524357 ,0.529011 ]
##maxP for generic system, t=1.0,tp=0.96,U=0.5,Up=0.48,complex_guass coef,B=0.1,localized in 6 sites
#maxp_Ulow_6=[1,0.959865,0.823407,0.760563,0.714737,0.683772,0.664105,0.634076,0.621026,0.621807,0.596454,0.572618,0.583776 ]

##S_EX with x-coarse of size 5 for generic system,t=1.0,tp=0.96,U=0.5,Up=0.48,complex_guass coef,B=0.1
L=np.array([5,10,15,20,25,30])
S0=[2.297,4.787,6.120,7.0388,7.7406603,8.30893817]
S=[1.905,3.239,  ,  ,  ]

##FOE with x-coarse of size 5
#L=np.array([5,10,15,20,25,30])
#foe=[  ,  , , , ,  ]

l=5
#l=6
Ne=3
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2=np.divide(N,np.power(M,2))

fig, ax = plt.subplots()
plt.xticks(np.arange(0, max(NM2), 5))
#ax.set_xticks(NM2)
#ax.set_xticklabels(case, minor=False,fontsize=8, rotation=90)

plt.scatter(NM2,maxp,s=30,color='g',marker="*",label='generic\ncomplex_gauss coef')
plt.scatter(NM2,maxp_U0,s=20,color='r',marker=".",label='U=Up=0\ncomplex_gauss coef')
#plt.scatter(NM2,maxp_real,s=15,color='b',marker="h",label='generic\nreal_gauss coef')
#plt.scatter(NM2,maxp_U0_real,s=15,color='m',marker=">",label='U=Up=0\nreal_gauss coef')
plt.scatter(NM2,maxp_Ulow_5,s=15,color='gray',marker=">",label='t=1.0,tp=0.96\nU=0.5,Up=0.48\ncomplex_gauss coef')
#plt.scatter(NM2,maxp_Ulow_6,s=15,color='m',marker=">",label='t=1.0,tp=0.96\nU=0.5,Up=0.48\ncomplex_gauss coef')


ax.set_axisbelow(True)
ax.minorticks_on()
# Customize the major grid
ax.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
# Customize the minor grid
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

# Turn off the display of all ticks.
ax.tick_params(which='both', # Options for both major and minor ticks
                top='off', # turn off top ticks
                left='on', # turn off left ticks
                right='off',  # turn off right ticks
                bottom='on') # turn off bottom ticks


plt.legend(loc='top right')
ax.text(34,0.7,'Beta=0.1\nM=10\n#particles=3')
#ax.text(7,0.85,'Beta=0.1\nM=20\n#particles=3')
plt.ylabel('maxP')
plt.xlabel('$\\frac{N}{M^{2}}$')
plt.show()

