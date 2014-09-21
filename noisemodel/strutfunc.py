from datatools.tempo import *
import numpy as np
from pylab import *
import matplotlib.pyplot as pyplot
fig = pyplot.figure()
ax1 = fig.add_subplot(2,1,1)

m = model('1713.Sep.dmx.par')
t = TOAfile('1713.Sep.T2.tim')
#m = model('1713.DMX.par')
#t = TOAfile('1713.8y.tim')
#m.thawall('DMX_')
#m.write('1713.DMX.par')
#sys.exit(0)
m.tempofit(t)
nm = model(m.newpar.parfile)
DMX, DMXErr, DMXR1, DMXR2 = nm.dmxlist
DMXR = {}
#print  DMXR1.keys()
#sys.exit(0)
for j in [k for k in DMXR1.keys() if not k == 1]:
    #DMXR[j] = (DMXR1[j] + DMXR2[j])/2
    toa_in_DMbin = []
    toa_sigma = []
    for i in range(len(m.toa)):
        if m.toa[i] > float(DMXR1[j]) and m.toa[i] < float(DMXR2[j]):
            toa_in_DMbin.append(m.toa[i])
            toa_sigma.append(m.err[i])
    if len(toa_in_DMbin) > 0:
        toa_in_DMbin = np.array(toa_in_DMbin)
        toa_sigma = np.array(toa_sigma)
        weights = 1/toa_sigma**2
        weighted_average = sum(weights*toa_in_DMbin)/sum(weights)
        DMXR[j] = weighted_average 
        #print len(toa_in_DMbin), weighted_average 
    else:
        DMXR[j] = (DMXR1[j] + DMXR2[j])/2

DMidx = []
for j in DMXR:
    DMidx.append((DMXR[j], j))
DMidx = np.array(DMidx, dtype=[('dm', float ), ('i', int)] )
DMidx.sort(order='dm')
idx = DMidx['i']

dmx = np.array([float(DMX[i]) for i in idx])
dmxerr = np.array([float(DMXErr[i]) for i in idx])
dmr = np.array([float(DMXR[i]) for i in idx])

ax1.errorbar(dmr, dmx, yerr=dmxerr, fmt='o')

Ntot =  len(dmx)
allpairs = np.array([(np.abs(dmr[j]-dmr[i]), (dmx[j]-dmx[i])**2, dmxerr[i]**2 + dmxerr[j]**2)  for i in range(Ntot-1) for j in range(i+1, Ntot)], dtype=[('delay', float), ('ddm', float), ('esq', float)])

allpairs.sort(order='delay')
tmax = allpairs['delay'].max()
tmin = allpairs['delay'].min()
bins = np.linspace(np.log10(tmin), np.log10(tmax+1),12)
#print bins
pairbin = []
i = 0
thisbin = []
for pair in allpairs:
    delay, ddmsq, errsq = pair
    if delay <= 10**bins[i+1]: 
        thisbin.append((delay, ddmsq, errsq))
    else:
        i+=1
        pairbin.append(np.array(thisbin, dtype=[('delay', float), ('ddm', float), ('esq', float)]))
        thisbin = [(delay, ddmsq, errsq)]
pairbin.append(np.array(thisbin, dtype=[('delay', float), ('ddm', float), ('esq', float)]))

#print pairbin[-1], 10**bins[-1]
binsize = []
for i in range(1,len(bins)):
    binsize.append(10**bins[i]-10**bins[i-1])
binsize = np.array(binsize)
#print pairbin
ax = fig.add_subplot(2,1,2)

delays = []
ddmsqs = []
ddmers = []
C = 4.148e3*1.e6
f = 1400.
fac = (2*np.pi*C/f)**2
for bin in pairbin:
    #ax.plot(bin['delay'],bin['ddm'] - bin['esq'], 'r.')
    delays.append(bin['delay'].mean())
    ddmsqs.append((bin['ddm'].mean() - bin['esq'].mean())*fac)
    #print bin['delay']
    #ddmsqs.append(bin['ddm'].mean()) #- bin[...,2].mean())
    #ddmers.append(bin['ddm'].std())
    #ddmers.append((2*np.sqrt(bin['ddm']*bin['esq'])).mean())
    ddmers.append((2*np.sqrt(sum(bin['ddm']*bin['esq'])))/bin['ddm'].size*fac)

#print ddmsqs
#print ddmers
#print len(delays), binsize.size
ax.errorbar(delays, ddmsqs, xerr=binsize/2, yerr=ddmers, fmt='o')
#ax.semilogx()
#ax.semilogy()
show()

data = np.vstack((delays, ddmsqs, binsize/2, ddmers))
np.save('strucfunc', data)
