from datatools.tempo import *
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.mlab as mlab
import scipy.stats

md = model('Feb.T2.RN.par')
tf = TOAfile('1713.Feb.T2.tim')
#print tf.groups.keys()
md.tempo2fit(tf)

padd_L = float(tf.toalist[tf.groups['ABPP-L'][0]].flags['padd'])
padd_S = float(tf.toalist[tf.groups['ABPP-S'][0]].flags['padd'])

"""
#print padd_L, padd_S
x = np.linspace(-0.0001, 0.0001, 30) + float(padd_L)

reschisq = []
for i in range(len(x)):
    for grp in ['ABPP-L', 'ABPP-S']:
        for k in tf.groups[grp]:
            toa = tf.toalist[k]
            toa.flags['padd'] = str(x[i])

    fo = open('1713.Mar.1d.tim', 'w')
    fo.write(tf.tempo2fmt())
    fo.close()

    del tf
    tf = TOAfile('1713.Mar.1d.tim')

    md.tempo2fit(tf)

    reschisq.append(md.chisq)

reschisq = np.array(reschisq)

np.savetxt('phasschisq.txt', (x, reschisq), fmt='%g')
#np.save('PJ2d', (X,Y,Z))
"""
#(X, Y, Z) = np.load('PJ2d.npy')

#plt.figure()
#plt.imshow(Z, origin='lower')
#plt.show()

x, y = np.loadtxt('phasschisq.txt')
ax = plt.subplot()
ax.plot(x, scipy.stats.distributions.chi2.cdf(y, md.dof), '-')
ax.set_xlabel('ABPP L&S phase jump')
ax.set_ylabel('$\chi^2 cdf$')
majorFormatter = FormatStrFormatter('%0.5f')
ax.xaxis.set_major_formatter(majorFormatter)
#majorFormatter = FormatStrFormatter('%5.0f')
#ax.yaxis.set_major_formatter(majorFormatter)
show()

#md.tempo2fit(tf)

#print tf.toalist[tf.groups['']]
#for toa in tf.toalist[]:

#print md.chisq, md.dof


#md.average()
#md.plot('mjd', 'averes')
#show()
