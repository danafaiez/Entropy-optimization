import matplotlib
import numpy as np;
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns; sns.set()
"""
#data_EVOLVMAX(Sxe):
data=[[0.211018],
[0.131963],
[0.205526],
[0.210418],
[0.122783],
[0.148105],
[0.144044],
[0.264051],
[0.125333],
[0.114635],
[0.144161],
[0.177961]]


#data_AVE(Sent):
data=[[0.167934],
[0.097548],
[0.203386],
[0.189947],
[0.161321],
[0.143618],
[0.135016],
[0.226040],
[0.199930],
[0.158678],
[0.198340],
[0.118242]]
"""

#data_OPTMAX(Sxe):
data=[[0.180303],
[0.157145],
[0.148879],
[0.161972],
[0.162782],
[0.143029],
[0.155761],
[0.216822],
[0.131606],
[0.152647],
[0.187246],
[0.201809]]


data = np.transpose(data)
fig = plt.figure()
subplots_adjust(left=None, bottom=None, right=None, top=0.45, wspace=None,hspace=0.2)
x=['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16']
label=['']
ax=sns.heatmap(data,cmap=plt.cm.get_cmap('PuBuGn'),cbar_kws={"orientation":"horizontal",'label': '<N>'},vmin=0, vmax=0.8)
ax.set_yticklabels(label, minor=False)
#ax.set_xlabel('Lattice site number',fontsize=10)
#ax.set_aspect('equal','box-forced')
ax.set_xticklabels(x)

#SxE heat map lines:
ax.vlines(np.arange(0,17,4), *ax.get_xlim(), color='dimgray',linewidth=4)
ax.vlines(np.arange(0,17,1), *ax.get_xlim(), color='dimgray',linewidth=1)
"""
#Sent heat map lines:
ax.vlines([0,4,16], *ax.get_xlim(), color='dimgray',linewidth=4)
ax.vlines(np.arange(0,16,1), *ax.get_xlim(), color='dimgray',linewidth=1)
"""
plt.savefig("heatxE_paper.svg",bbox_inches='tight')
plt.show()
