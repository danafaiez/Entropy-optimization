import matplotlib
import numpy as np;
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns; sns.set()
"""
#data_aveSent:
data=[[0.094224],
[0.141048],
[0.211777],
[0.172796],
[0.144082],
[0.090486],
[0.052831],
[0.039195],
[0.077197],
[0.158195],
[0.232840],
[0.217493],
[0.174628],
[0.100477],
[0.061765],
[0.030964]]
"""

#data_maxSent:
data=[[0.085437],
[0.113153],
[0.148205],
[0.193699],
[0.205212],
[0.223058],
[0.178883],
[0.120615],
[0.086461],
[0.090028],
[0.106084],
[0.131608],
[0.120486],
[0.104635],
[0.059090],
[0.033347]] 




data = np.transpose(data)
fig = plt.figure()
#subplots_adjust(left=None, bottom=None, right=None, top=0.45, wspace=None,hspace=0.2)
x=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16']
label=['']
ax=sns.heatmap(data,cmap=plt.cm.get_cmap('PuBuGn'),cbar_kws={"orientation":"horizontal",'label': '<N>'},vmin=0, vmax=0.8)
ax.set_yticklabels(label, minor=False)
#ax.set_xlabel('Lattice site number',fontsize=10)
#ax.set_aspect('equal','box-forced')
ax.set_xticklabels(x)

#Sent heat map lines:
ax.vlines([0,4,16], *ax.get_xlim(), color='dimgray',linewidth=4)
ax.vlines(np.arange(0,16,1), *ax.get_xlim(), color='dimgray',linewidth=1)
plt.savefig("heat_ave_max_Sent.svg",bbox_inches='tight')
plt.show()
