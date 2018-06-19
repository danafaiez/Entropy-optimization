#plotting max prob of localization in the 1st 5 sites,for a generic case, with beta=0, Ne=3, for different L
#!usr/local.bin/python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
#from matplotlib.offsetbox import AnchoredText
L=[10,12,15,17,20,23,25,27,30]
max=[0.73,0.68,0.65,0.635,0.59,0.584,0.563,0.567,0.54]
fig, ax = plt.subplots()
ax.set_xticks(L)
#ax.set_xticklabels(L, minor=False)

plt.scatter(L,max,s=20, color='m', marker="*")
plt.ylabel('maxP')
plt.xlabel('L')
#plt.legend(loc='lower left')
#plt.text(0, 1.51,r'integrable \n g-:#bath_sites=8')
plt.show()

