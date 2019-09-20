import matplotlib
import numpy as np;
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns; sns.set()

#data_Sxe:
data=[[0.032140],
[0.030747],
[0.036591],
[0.028444],
[0.403696],
[0.437729],
[0.416554],
[0.394815],
[0.033245],
[0.040269],
[0.025550],
[0.018386],
[0.029891],
[0.019200],
[0.027118],
[0.025624]]
"""
#data_Sent:
data=[[0.000224],
[0.000211],
[0.000225],
[0.000126],
[0.175656],
[0.202968],
[0.193706],
[0.181722],
[0.101337],
[0.212285],
[0.165999],
[0.147534],
[0.210366],
[0.193470],
[0.088930],
[0.125242]]
"""
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
