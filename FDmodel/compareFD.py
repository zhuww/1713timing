from pylab import *
from datatools.tempo import *
import numpy as np
secperday = 24*3600

ax1 = subplot(211)
ax2 = subplot(212)

m = model('1713.FD0.par')
#m = model('1713.pbd.par')
t = TOAfile('1713.8y.tim')
m.tempofit(t)
#m.freezeall('DM')
#m.write('1713.DM0.par')
m.groups['GUPPI'] = m.groups['GUPPI-L'] + m.groups['GUPPI-8']
m.groups['PUPPI'] = m.groups['PUPPI-L'] + m.groups['PUPPI-S']
m.average(freqbin=20.)

colors = {'GUPPI':'r', 'PUPPI':'c'}
m.plot('freq','averes', groups=['GUPPI', 'PUPPI'], colors={'GUPPI':'r', 'PUPPI':'c'}, LegendOn=True, LegendLoc=1, capsize=0, elinewidth=1, marker='.', ms=0, NoZeroLine=True, ax=ax1)

m1=model('1713.DM0.par')
FD1 = float(m1.FD1[0])
FD2 = float(m1.FD2[0])
FD3 = float(m1.FD3[0])
FD4 = float(m1.FD4[0])
p = (FD1, FD2, FD3, FD4)

def func(f, C):
    logf = np.log(f/1000.)
    return (p[0] * logf + p[1]*logf**2 + p[2]*logf**3 + p[3]*logf**4)*1.e6 + C

import scipy.optimize as opt

popt, pcov = opt.curve_fit(func, m.freq, m.res)
C = popt[0]

freq = np.linspace(700, 2500, 50)
freq.sort()
ax1.plot(freq, func(freq, C),'-')

del m,t

m = model('1713.pbd.par')
#m = model('1713.pbd.par')
t = TOAfile('1713.8y.tim')
m.tempofit(t)
m.groups['GUPPI'] = m.groups['GUPPI-L'] + m.groups['GUPPI-8']
m.groups['PUPPI'] = m.groups['PUPPI-L'] + m.groups['PUPPI-S']
m.average(freqbin=20.)
colors = {'GUPPI':'r', 'PUPPI':'c'}
m.plot('freq','averes', groups=['GUPPI', 'PUPPI'], colors={'GUPPI':'r', 'PUPPI':'c'}, LegendOn=True, LegendLoc=1, capsize=0, elinewidth=1, marker='.', ms=0, NoZeroLine=True, ax=ax2)
ax2.axhline(y=0, linestyle='--')

show()
