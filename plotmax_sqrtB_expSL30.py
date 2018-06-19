#L=30, Ne=3, generic, using newPsi unitary method to maximize the probability of localization in the first 5 or 10 or 6 sites vs sqrt(Beta) or exp(Delta(S)).
#!usr/local.bin/python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#generic case, localizing in the first 5 sites
#Beta=[0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.0]
Beta=[0,0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0]
sqrt_Beta=np.sqrt(Beta)
maxp=[0.544634,0.448685 ,0.288223,0.138452 ,0.053152,0.017839,0.005653,0.001780,0.000575,0.000194,0.000069 ]
#Beta_used_generic=[3.6998e-04, ,2.00286742e-1 ,4.0069205e-1 ,6.004501e-1,8.0778565e-1 ,1.036440 ,1.219949130 ,1.3791 ,1.758184 , 1.758184 ,1.758184]
S=np.array([8.3089375328 ,8.06803376 ,7.3239987,6.280596114 ,5.21391032781 ,4.282244123 ,3.73269715414 , 3.36635776 ,2.753692,2.7536926,2.7536926])
expDeltaS=np.exp([S-S[0]])
expDeltaS=list(expDeltaS)

##S_EX with x-coarse of size 5
#S=[  , ,  ,  ,  ,  ]

##FOE with x-coarse of size 5
##foe=[  ,  , , , ,  ]

##generic, t=U=1.0,tp=Up=0.96 localizing in the first 10 sites (M=120)  
maxp_10=[0.723383,0.700778,0.633898,0.520311,0.388789,0.270550,0.179787,0.117113,0.076026,0.049653,0.032737]

##generic, t=1.0,tp=0.96 U=0.5,Up=0.48 localizing in the first 6 sites  
maxp_Ulow_6=[0.598496,0.547969,0.412533,0.245610,0.117935,0.049086,0.019120,0.007365,0.002900,0.001188,0.00051]

##generic, t=U=1.0,tp=Up=0.96 localizing in the first 5 sites, real-gaussian coef
maxp_realG_coef=[0.442079,0.355797,0.221225,0.104470,0.039826,0.013409,0.004263,0.001339,0.000428 ,0.000142,0.000050]


##generic,t=1.0,tp=0.96 U=0.5,Up=0.48 localizing in the first 5 sites  
maxp_Ulow_5=[0.553522,0.482010,0.329270,0.168357,0.067853,0.023646,0.007726,0.002497,0.000823,0.000282,0.000102]

##generic,t=1.0,tp=0.96 U=Up=0, localizing in the first 5 sites  
maxp_U0=[0.538325,0.491957,0.351423,0.193360,0.087653,0.033710,0.012117,0.004283]

##generic, t=1.3,tp=0.96,U=0.50,Up=0.48 localizing in the first 5 sites
maxp_lU_ht=[0.547476,0.459427,0.268296,0.108295,0.034217,0.009773,0.002769,0.000816,0.000256,0.000086,0.000031]

##generic, t=1.3,tp=1.30,U=0.50,Up=0.48 localizing in the middle 5 sites
maxp_lU_http_m5=[0.582759,0.486882,0.280570,0.108982,0.035328,0.011419,0.004038,0.001622, 0.000745,0.000389,0.000227]

##generic, t=1.9,tp=1.90,U=0.50,Up=0.48 localizing in the middle 5 sites
maxp_lU_hhttp_m5=[0.592194,0.409526,0.126771,0.025407, 0.005394,0.001460,0.000520,0.000236,0.000129, 0.000081,0.000056] 

##generic, t=1.3,tp=1.30,U=0.50,Up=0.48 localizing in the end 5 sites
maxp_lU_http_e5=[0.555551,0.442409,0.223202,0.072557,0.018237,0.004219,0.000992,0.000249,0.000068,0.000020,0.000007]

##Integrable, t=1,tp=0.00,U=0.5,Up=0 localizing in the first 6 sites
maxp_int_lowU_6=[0.590269,0.568529,0.528301,0.468844,0.406139,0.336265,0.271645,0.215479,0.168844,0.131307,0.101694]

##Integrable, t=1,tp=0.00,U=0.5,Up=0 localizing in the first 5 sites
maxp_int_lowU_5=[0.540194,0.512869,0.455390,0.378262,0.296576,0.222333,0.161367,0.114582,0.080311,0.055930 ,0.038879]

##Integrable, t=1.3,tp=0.00,U=0.5,Up=0 localizing in the first 5 sites
maxp_int_lU_hK_5=[0.541545,0.499928,0.420566,0.313175,0.216564,0.142778,0.091613,0.058060,0.036694,0.034217,0.014846]

##Integrable, t=1.9,tp=0.00,U=0.5,Up=0 localizing in the first 5 sites
maxp_int_lU_hhK_e5=[0.553906,0.478838,0.332228,0.196595,0.106571,0.055921,0.029228,0.015440,0.008300,0.004556,0.002558]

##Integrable, t=1.9,tp=0.00,U=0.5,Up=0 localizing in the middle 5 sites
maxp_int_lU_hhK_m5=[0.572635,0.506867,0.374892,0.245721,0.152472,0.094121,0.059328,0.038559,0.025878,0.017899,0.012717]

##Pmax for t=1,U=0.50,Up=0.48 system at B=0, localizing in the first 5 sites, vs tp
tp=[0,0.2,0.48,1,1.3,1.5]
maxp_vs_tp=[0.555505,0.565333,0.560498,0.552996,0.555709,0.553151] 

##Pmax for tp=1,U=0.50,Up=0.48 system at B=0, localizing in the first 5 sites, vs t
t=[0.5,1,1.3,1.5]
maxp_vs_t=[0.558871,0.552996,0.542269,0.528590]

##Beta_used=[1.628444e-04,2.0071e-01,4.002770,6.00666,7.97378,1.0216,1.272872,1.43995]
S_U0=np.array([8.3089381,8.0775766,7.3798250,6.391861,5.4094662,4.481239,3.729836,3.36336])
expDeltaS_U0=np.exp([(S_U0)-(S_U0[0])])
expDeltaS_U0=list(expDeltaS_U0)



fig, ax = plt.subplots()
##plotting Pmax vs tp,B=0 localized in 5 sites
#plt.scatter(tp,maxp_vs_tp,s=25,color='g', marker="<",label='L=30\nsize of box=5\n#particles=3\nB=0\nt=1.0\nU=0.5,Up=0.48\ncomplexgaussian coef')
##plotting Pmax vs t,B=0 localized in 5 sites
#plt.scatter(t,maxp_vs_t,s=25,color='r', marker="<",label='L=30\nsize of box=5\n#particles=3\nB=0\ntp=1.0\nU=0.5,Up=0.48\ncomplexgaussian coef')
##plotting the generic case - localized in 5 sites
#plt.scatter(sqrt_Beta,maxp,s=30,color='m', marker="o",label='L=30\nsize of box=5\n#particles=3\nGeneric system:\nt=1.0,tp=0.96\nU=1.0,Up=0.96\ncomplex gaussian coef')
#plt.scatter(expDeltaS,maxp,s=30,color='g', marker="o",label='L=30\nsize of box=5\n#particles=3')
##plotting the U=Up=0 case - localized in 5 sites
#plt.scatter(expDeltaS_U0,maxp_U0,s=30,color='b', marker="o",label='U=Up=0\nL=30\nsize of box=5\n#particles=3')
##plotting the generic case - localized in 10 sites
#plt.scatter(expDeltaS,maxp_10,s=30,color='gray', marker="o",label='L=30\nsize of box=10\n#particles=3')
##plotting the generic case - localized in 5 sites - real gaussian coef
#plt.scatter(sqrt_Beta,maxp_realG_coef,s=30,color='g', marker="o",label='L=30\nsize of box=5\n#particles=3\nreal gaussian coef')
##plotting the generic case - U=0.5,Up=0.48 - localized in 5 sites - complex gaussian coef
#plt.scatter(sqrt_Beta,maxp_Ulow_5,s=25,color='m',marker="o",label='L=30\nsize of box=5\n#particles=3\nGeneric system:\nt=1.0,tp=0.96\nU=0.5,Up=0.48\ncomplex gaussian coef')
##plotting the integrable case - t=1.0,tp=0.00,U=0.5,Up=0.00 - localized in 5 sites - complex gaussian coef
#plt.scatter(sqrt_Beta,maxp_int_lowU_5,s=30,color='r',marker="o",label='L=30\nsize of box=5\n#particles=3\nIntegrable system:\nt=1.0,tp=0.00\nU=0.5,Up=0.00\ncomplex gaussian coef')
##plotting the generic case - U=0.5,Up=0.48 - localized in 6 sites - complex gaussian coef
#plt.scatter(sqrt_Beta,maxp_Ulow_6,s=25,color='b',marker="o",label='L=30\nsize of box=6\n#particles=3\nGeneric system:\nt=1.0,tp=0.96\nU=0.5,Up=0.48\ncomplex gaussian coef')
##plotting the integrable case - t=1.0,tp=0.00,U=0.5,Up=0.00 - localized in 6 sites - complex gaussian coef
#plt.scatter(sqrt_Beta,maxp_int_lowU_6,s=30,color='r',marker="o",label='L=30\nsize of box=6\n#particles=3\nIntegrable system:\nt=1.0,tp=0.00\nU=0.5,Up=0.00\ncomplex gaussian coef')
##plotting Pmax vs sqrtB, localization on first 5 sites, t=1.3,tp=1.30,U=0.50,Up=0.48
#plt.scatter(sqrt_Beta,maxp_lU_http_m5,s=25,color='m',marker="<",label='L=30\nlocalizing in middle 5 sites\n#particles=3\nt=1.3,tp=1.30\nU=0.5,Up=0.48\ncomplex gaussian coef')
plt.scatter(sqrt_Beta,maxp_lU_http_m5,s=25,color='m',marker="<",label='loc:middle 5 sites\nt=1.3,tp=1.30\nU=0.5,Up=0.48')

##plotting Pmax vs sqrtB, localization on middle 5 sites,t=1.3,tp=1.30,U=0.50,Up=0.48 
#plt.scatter(sqrt_Beta,maxp_lU_http_e5,s=25,color='b',marker="<",label='L=30\nlocalizing in first 5 sites\n#particles=3\nt=1.3,tp=1.30\nU=0.5,Up=0.48\ncomplex gaussian coef')
plt.scatter(sqrt_Beta,maxp_lU_http_e5,s=25,color='b',marker="<",label='loc:first 5 sites\nt=1.3,tp=1.30\nU=0.5,Up=0.48')

##plotting Pmax vs sqrtB, generic, t=1.9,tp=1.90,U=0.50,Up=0.48 localizing in the middle 5 sites
#plt.scatter(sqrt_Beta,maxp_lU_hhttp_m5,s=25,color='r',marker="<",label='L=30\nlocalizing in middle 5 sites\n#particles=3\nt=1.9,tp=1.90\nU=0.5,Up=0.48\ncomplex gaussian coef')
plt.scatter(sqrt_Beta,maxp_lU_hhttp_m5,s=25,color='r',marker="<",label='loc.:middle 5 sites\nt=1.9,tp=1.90\nU=0.5,Up=0.48')

##plotting Pmax vs sqrtB,integrable, t=1.9,tp=1.90,U=0.50,Up=0.48 localizing in the first 5 sites
#plt.scatter(sqrt_Beta,maxp_int_lU_hhK_e5,s=25,color='g',marker="o",label='L=30\nlocalizing in first 5 sites\n#particles=3\nt=1.9,tp=0.00\nU=0.5,Up=0.00\ncomplex gaussian coef')
plt.scatter(sqrt_Beta,maxp_int_lU_hhK_e5,s=25,color='g',marker="o",label='loc.:first 5 sites\nt=1.9,tp=0.00\nU=0.5,Up=0.00')

##plotting Pmax vs sqrtB,integrable, t=1.9,tp=1.90,U=0.50,Up=0.48 localizing in the middle 5 sites
#plt.scatter(sqrt_Beta,maxp_int_lU_hhK_m5,s=25,color='gray',marker="o",label='L=30\nlocalizing in middle 5 sites\n#particles=3\nt=1.9,tp=0.00\nU=0.5,Up=0.00\ncomplex gaussian coef')
plt.scatter(sqrt_Beta,maxp_int_lU_hhK_m5,s=25,color='gray',marker="o",label='loc.:middle 5 sites\nt=1.9,tp=0.00\nU=0.5,Up=0.00')

#plt.xticks(np.arange(0, 2, 0.1))

plt.yticks(np.arange(0,0.75,0.1))
plt.legend(loc='bottom left')
#plt.legend(loc='top right')
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

#plt.xlabel('expDeltaS')
plt.xlabel('$\\sqrt{\\beta}$')
#plt.ylabel('$\\frac{maxP(Beta)}{maxP(Beta=0)}$')
plt.ylabel('Pmax')
plt.show()

