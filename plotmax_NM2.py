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
maxp_Ulow_6=[1,0.959865,0.823407,0.760563,0.714737,0.683772,0.664105,0.634076,0.621026,0.621807,0.596454,0.572618,0.583776 ]

##maxP for generic system, t=1.9,tp=1.90,U=0.5,Up=0.48,complex_guass coef,B=0.01,localized in middle 5 sites
L5=np.array([5,10,12,15,18,20,22,25,28,30])
maxp_Ulow_m5=[1,0.745414,0.720627,0.616848, 0.611636,0.613332,0.602319 ,0.545235,0.589383 ,0.589317]
maxp_Ulow_e5=[1,0.717694,0.694804,0.642354,0.602576,0.601286,0.571663,0.568835,0.560111,0.545582]
##S_ex corresponding to the state with Pmax at t=1.9,tp=1.90,U=0.5,Up=0.48,complex_guass coef,B=0.01,localized in middle 5 sites 
L5S=np.array([5,10,15,20,25,30])
#???S=[1.349182,4.020062,3.774280,5.244437,4.714170,5.916198]

##Average S over time (t=2000,delta_t=4)
#Save=[2.021148,4.353730105,5.701217972,6.615725676,7.3150042484,7.88537424]
#R=np.divide(Save,S)
##FOE with x-coarse of size 5
#L=np.array([5,10,15,20,25,30])
#foe=[  ,  , , , ,  ]


##maxP for generic system, t=1.9,tp=1.90,U=0.5,Up=0.48,complex_guass coef,B=0.01,localized in middle 6 sites
maxp_Ulow_m6=[1,0.796813,0.769364,0.720122,0.642048,0.639189,0.624886,0.639288,0.609841,0.605166]

##maxP for generic system, t=1.9,tp=1.90,U=0.5,Up=0.48,complex_guass coef, B=0.1,localized in middle 5 sites coef
#L5=np.array([5,10,12,15,18,20,22,25,28,30])
#maxp_Ulow_m5=[1 ,0.746056 ,0.701287 ,0.565115,0.567391 ,0.570588 ,0.555708 ,0.484568 ,0.539691,0.535790]
#SL5=[1.192989,3.950510,4.540390 ,4.795842,5.067436 ,0.555708 ,4.837910 ,3.320505 ,3.314089 ]



l=5
#l=6
Ne=3
N=factorial(L)/(factorial(Ne)*factorial(L-Ne))
M=factorial(l)/(factorial(Ne)*factorial(l-Ne))
NM2=np.divide(N,np.power(M,2))

N_L5=factorial(L5)/(factorial(Ne)*factorial(L5-Ne))
NM2_L5=np.divide(N_L5,np.power(M,2))
fig, ax = plt.subplots()
#plt.xticks(np.arange(0, max(NM2), 5))

plt.xticks(np.arange(0, max(NM2_L5), 5))
#ax.set_xticks(NM2)
#ax.set_xticklabels(case, minor=False,fontsize=8, rotation=90)

#plt.scatter(NM2,maxp,s=30,color='g',marker="*",label='generic\ncomplex_gauss coef')
#plt.scatter(NM2,maxp_U0,s=20,color='r',marker=".",label='U=Up=0\ncomplex_gauss coef')
#plt.scatter(NM2,maxp_real,s=15,color='b',marker="h",label='generic\nreal_gauss coef')
#plt.scatter(NM2,maxp_U0_real,s=15,color='m',marker=">",label='U=Up=0\nreal_gauss coef')
#plt.scatter(NM2,maxp_Ulow_5,s=15,color='gray',marker=">",label='t=1.0,tp=0.96\nU=0.5,Up=0.48\ncomplex_gauss coef')
#plt.scatter(NM2,maxp_Ulow_6,s=15,color='m',marker=">",label='t=1.0,tp=0.96\nU=0.5,Up=0.48\ncomplex_gauss coef')
#plt.scatter(NM2_L5,maxp_Ulow_m5,s=15,color='g',marker=">",label='t=1.9,tp=1.9\nU=0.5,Up=0.48\ncomplex_gauss coef\nMiddle;M=10')
plt.scatter(NM2_L5,maxp_Ulow_m5,s=15,color='g',marker=">",label='middle 5 sites')


#plt.scatter(NM2_L5,maxp_Ulow_m6,s=15,color='b',marker=">",label='t=1.9,tp=1.9\nU=0.5,Up=0.48\ncomplex_gauss coef\nMiddle;M=20')
plt.scatter(NM2_L5,maxp_Ulow_m6,s=15,color='b',marker=">",label='middle 6 sites')

#plt.scatter(NM2_L5,maxp_Ulow_e5,s=15,color='r',marker=">",label='t=1.9,tp=1.9\nU=0.5,Up=0.48\ncomplex_gauss coef\nSide;M=10')
plt.scatter(NM2_L5,maxp_Ulow_e5,s=15,color='r',marker=">",label='first 5 sites')
#plt.scatter(L5S,R,s=15,color='r',marker="m",label='t=1.9,tp=1.9\nU=0.5,Up=0.48\ncomplex_gauss coef\nM10;size_box=5')
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
ax.text(34,0.92,'Beta=0.01\n#particles=3')
#ax.text(7,0.85,'Beta=0.1\nM=20\n#particles=3')
plt.ylabel('maxP')
plt.xlabel('$\\frac{N}{M^{2}}$')
plt.show()

