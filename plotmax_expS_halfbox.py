#L=20, Ne=4, generic, using newPsi unitary method to maximize the probability of localization in the first 10 sites vs exp{-[S(beta)-S(0)]}.
#!usr/local.bin/python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import *
b=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]

S_beta=np.array([8.485701086,8.417177326,8.2010516626,7.83327707,7.3293718668,
6.73324186526,6.075724,5.41890,4.7730014,4.14346332,3.654101227,3.148837232 ,2.78569079,
2.321609774,2.044814587 ,1.7217926,1.7217926,1.105906,1.105906,1.105906,1.105906])
expDeltaS=np.exp([S_beta-S_beta[0]])
expDeltaS=list(expDeltaS)
S_beta=list(S_beta)
maxp=np.array([0.763794,0.747246,0.726650,0.683923,0.620971,0.536772,0.441354,0.343988,
0.254306,0.180415,0.124234,0.083540,0.055473,0.036655,0.024248 ,0.016129,0.010823,0.007343,
0.005046,0.003515,0.002484])
sqrt_beta=np.sqrt(b)

maxp_ratio=maxp/0.763794
maxp=list(maxp)
maxp_ratio=list(maxp_ratio)
fig, ax = plt.subplots()
#ax.set_xticks(expDeltaS)


#plt.scatter(expDeltaS,maxp,s=30,color='r', marker="o",label='L=20\nsize of box=10\n#particles=4')
#plt.scatter(expDeltaS,maxp_ratio,s=30,color='b', marker="o",label='L=20\nsize of box=10\n#particles=4')
plt.scatter(sqrt_beta,expDeltaS,s=30,color='r', marker="o",label='L=20\nsize of box=10\n#particles=4')
#plt.scatter(b,S_beta,s=30,color='r', marker="o",label='L=20\nsize of box=10\n#particles=4')

plt.legend(loc='top right')
#plt.ylabel('maxP(Beta)')
#plt.ylabel('$\\frac{maxP(Beta)}{maxP(Beta=0)}$')
plt.xlabel('$\\sqrt{\\beta}$')
plt.ylabel('expDeltaS')
#plt.xlabel('expDeltaS')
plt.show()

