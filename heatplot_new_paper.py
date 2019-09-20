import matplotlib
import numpy as np;
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pylab import *
from scipy import stats
from scipy.optimize import curve_fit
import seaborn as sns; sns.set()

#data_Sxe:
data=[[0.026942],
[0.029773],
[0.038699],
[0.052995],
[0.339777],
[0.342788],
[0.314807],
[0.315846],
[0.134389],
[0.068634],
[0.131681],
[0.076173],
[0.033919],
[0.019732],
[0.048143],
[0.025700]]

"""
#data_Sent:
data=[[0.000516],
[0.000372],
[0.000662],
[0.000552],
[0.112562],
[0.139109],
[0.176077],
[0.166501],
[0.209886],
[0.150686],
[0.195781],
[0.139608],
[0.204502],
[0.227058],
[0.158503],
[0.117625]]
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
