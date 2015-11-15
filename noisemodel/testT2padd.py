from datatools.tempo import *
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import FormatStrFormatter
import matplotlib.mlab as mlab

md = model('Feb.T2.RN.par')
tf = TOAfile('1713.Feb.T2.tim')
#print tf.groups.keys()

padd_L = float(tf.toalist[tf.groups['ABPP-L'][0]].flags['padd'])
padd_S = float(tf.toalist[tf.groups['ABPP-S'][0]].flags['padd'])

#"""
#print padd_L, padd_S
#errors = np.linspace(-0.00003, 0.00003, 30)
x = np.linspace(-0.00003, 0.00003, 10) + padd_L
y = np.linspace(-0.0001, 0.0001, 20) + padd_S
X, Y = np.meshgrid(x,y)
Z = np.zeros((len(x), len(y)))

reschisq = []
for i in range(len(x)):
    for k in tf.groups['ABPP-L']:
        toa = tf.toalist[k]
        toa.flags['padd'] = str(x[i])
    for j in range(len(y)):
        for l in tf.groups['ABPP-S']:
            toa = tf.toalist[l]
            toa.flags['padd'] = str(y[j])

        fo = open('1713.Mar.test.tim', 'w')
        fo.write(tf.tempo2fmt())
        fo.close()

        del tf
        tf = TOAfile('1713.Mar.test.tim')

        md.tempo2fit(tf)

        #reschisq.append(md.chisq)
        Z[i,j] = md.chisq

#reschisq = np.array(reschisq)

#np.savetxt('phasschisq.txt', (errors, reschisq), fmt='%g')
np.save('PJ2d', (X,Y,Z))
#"""
#(X, Y, Z) = np.load('PJ2d.npy')

plt.figure()
#CS = plt.contour(X, Y, Z)
plt.imshow(Z, origin='lower')
plt.show()

#x, y = np.loadtxt('phasschisq.txt')
#ax = plt.subplot()
#ax.plot(float(padd_L) + x, y, '-')
#ax.set_xlabel('ABPP L&S phase jump')
#ax.set_ylabel('$\chi^2$')
#majorFormatter = FormatStrFormatter('%0.5f')
#ax.xaxis.set_major_formatter(majorFormatter)
#majorFormatter = FormatStrFormatter('%5.0f')
#ax.yaxis.set_major_formatter(majorFormatter)
#show()

#md.tempo2fit(tf)

#print tf.toalist[tf.groups['']]
#for toa in tf.toalist[]:

#print md.chisq, md.dof


#md.average()
#md.plot('mjd', 'averes')
#show()
