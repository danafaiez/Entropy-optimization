#L=20, Ne=2, generic, using newPsi unitary method to maximize the probability of localization in the first x sites at different beta.
#!usr/local.bin/python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

beta=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
maxp_loc3=[0.512804,0.450350,0.316851,0.176871,0.084650,0.038262,0.017494,0.008386,0.004274,0.002320,0.001338]
maxp_loc4=[0.635786,0.580528,0.447605,0.305097,0.179309,0.099716,0.055189,0.031345,0.018507,0.011405,0.007332]
maxp_loc5=[0.700648,0.668773,0.574928,0.433370,0.293794,0.188852,0.120098,0.077378,0.051109,0.034758,0.024353]

sqrt_beta=np.sqrt(beta)
fig, ax = plt.subplots()
ax.set_xticks(beta)
#ax.set_xticklabels(case, minor=False,fontsize=8, rotation=90)
plt.scatter(sqrt_beta,maxp_loc3,s=30,color='r', marker="*",label='maximum probability of localization in the first 3 sites')
plt.scatter(sqrt_beta,maxp_loc4,s=30,color='m', marker="h",label='maximum probability of localization in the first 4 sites')
plt.scatter(sqrt_beta,maxp_loc5,s=30,color='b', marker="+",label='maximum probability of localization in the first 5 sites')
plt.legend(loc='bottom left')
#ax.text(1.05,0.55, r'L=20 #particles=2', fontsize=8)
plt.ylabel('P$_{\mathrm{max}}$', fontsize=12)
plt.xlabel('$\sqrt{\\beta}$')
plt.savefig("PmaxvsBeta.eps",format='eps',bbox_inches='tight')
plt.show()

