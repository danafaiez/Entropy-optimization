#generic, using newPsi unitary method to maximize the probability of localization in the first L/2 .
#plot pf maxP vs L, for B=0.5 and B=2. With Ne=3.
#!usr/local.bin/python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

L=[10,12,14,16,18,20,22,24,26,28,30]
maxp_beta_half=[0.405210,0.514098,0.581179,0.593800,0.662096,0.700631,0.699967,0.748498,0.765690,0.782113,0.791172 ]
maxp_beta2=[0.006651,0.007894,0.013228,0.020264,0.033395,0.054355,0.082158,0.119660,0.161418,0.209949,0.264791]

fig, ax = plt.subplots()
ax.set_xticks(L)
#ax.set_xticklabels(case, minor=False,fontsize=8, rotation=90)
plt.scatter(L,maxp_beta_half,s=30,color='r', marker="o",label='Beta=0.5')
plt.scatter(L,maxp_beta2,s=30,color='m', marker="h",label='Beta=2')
#plt.scatter(beta,maxp_loc5,s=30,color='b', marker="+",label='maximum probability of localization in the first 5 sites')
plt.legend(loc='top right')
plt.ylabel('maxP')
plt.xlabel('L')
ax.text(9.35,0.65, r'#particles=3 , box=L/2', fontsize=8)
plt.show()

