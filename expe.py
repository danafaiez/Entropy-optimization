#!usr/local.bin/python
import matplotlib
#matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
from pylab import *

#non-integrable t=tp=1.9,U=Up=0.5, L=30, Ne=2, B=0.1
bath=[3,4,5,7,9,11,13,15,17,30]
#Ei=[-0.242542,-1.231589,-1.445467,-1.738883,-2.033278,-2.062927,-2.422980,-2.521376,-2.698577,-3.038908]
#Eo=[-2.207473,-1.827231,-1.520376,-1.131927,-0.993362,-0.813000,-0.536509,-0.444356,-0.275174,0]
#maxEi=[-0.242542,-0.932962,-1.445467,-1.738883,-2.033278,-2.062927,-2.422980,-2.521376,-2.698577,-3.038908]
#Etot=list( map(add, Ei, Eo))


#non-integrable t=tp=1.9,U=Up=0.5, L=30, Ne=3, B=0.1
Ei=[ ]
Eo=[ ]
Etot=list( map(add, Ei, Eo))


fig, ax = plt.subplots()
plt.scatter(bath,Ei,s=25,color='m',marker="o",label='Ei')
plt.scatter(bath,Eo,s=25,color='b',marker="*",label='Eo')
plt.scatter(bath,Etot,s=25,color='y',marker=".",label='Etot')
#plt.scatter(bath,maxEi,s=20,color='black',marker=".",label='maxEi')
#ax.text(21,-0.5,'L=30\n#particles=2\nL=30,t=tp=1.9,U=Up=0.5,B=0.1')
ax.text(21,-0.5,'L=30\n#particles=3\nL=30,t=tp=1.9,U=Up=0.5,B=0.1')

plt.legend(loc='top right')
show()
