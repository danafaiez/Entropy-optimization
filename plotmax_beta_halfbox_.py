#L=20, Ne=4, generic, using newPsi unitary method to maximize the probability of localization in the first 10 sites vs sqrt(beta).
#!usr/local.bin/python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
#maxp_beta_half=[0.405210,0.514098,0.581179,0.593800,0.662096,0.700631,0.699967,0.748498,0.765690,0.782113,0.791172 ]
#maxp_beta2=[0.006651,0.007894,0.013228,0.020264,0.033395,0.054355,0.082158,0.119660,0.161418,0.209949,0.264791]


beta=[0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
sqrt_beta=np.sqrt(beta)
maxp=[0.763794,0.747246,0.726650,0.683923,0.620971,0.536772,0.441354,0.343988,0.254306,0.180415,0.124234,0.083540,
0.055473,0.036655,0.024248,0.016129,0.010823,0.007343,0.005046,0.003515,0.002484]
fig, ax = plt.subplots()
plt.xticks(np.arange(0, 2, 0.1))
plt.yticks(np.arange(0, 1, 0.1))
plt.scatter(sqrt_beta,maxp,s=30,color='b', marker="h",label='L=20\n#particles=4\nsize of box=10')

#plt.scatter(beta,maxp,s=30,color='m',marker="h")

plt.legend(loc='top right')
plt.ylabel('maxP')
plt.xlabel('$\sqrt{\\beta}$')

ax.set_axisbelow(True)
ax.minorticks_on()
# Customize the major grid
ax.grid(which='major', linestyle=':', linewidth='0.5', color='gray')
# Customize the minor grid
ax.grid(which='minor', linestyle=':', linewidth='0.5', color='black')

# Turn off the display of all ticks.
ax.tick_params(which='both', # Options for both major and minor ticks
                top='off', # turn off top ticks
                left='off', # turn off left ticks
                right='off',  # turn off right ticks
                bottom='on') # turn off bottom ticks

plt.show()

