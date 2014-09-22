from datatools.tempo import *
from pylab import *
import numpy as np
from math import *
year = 365.24218967

m = model('1713.Sep.KOM.par')
#t = TOAfile('1713.Nov.Paul.tim')
#m.tempofit(t)
#m.average(groups='allinone')

#wrms = m.wrms['all']
wrms = 0.212

CPMSP_sc10= lambda nv, nvdot, T:exp(1.6) * nv**-1.4*abs(nvdot/1.e-15)**1.1*T**2

#Hobbs = lambda nv, nvdot, T:10.**-11.5 * nv*-0.4 * abs(nvdot/1.e-15)**0.8
Hobbs = lambda nv, nvdot, T:10**(-1.37*log10(nv**0.29*abs(nvdot)**-0.55) + 0.52)

Tspan = float(m.FINISH - m.START)/year
nv = float(m.F0[0])
nvdot = float(m.F1[0])

print wrms
print CPMSP_sc10(nv, nvdot, Tspan)
print Hobbs(nv, nvdot, Tspan)

#""" 
#The timing noise calculations:
#"""

#groups = m.averes.keys()
#groups.sort()
#for grp in groups:
    #errorbar(m.avetoa[grp], m.averes[grp], yerr=m.aveerr[grp], fmt = '.', mew=0, label=grp)
#xlabel(r'MJD (day)')
#ylabel(r'Residual ($\mu$s)')
#show()

#alltoa = []
#allres = []
#allerr = []
#for grp in groups:
    #alltoa.extend(m.avetoa[grp])
    #allres.extend(m.averes[grp])
    #allerr.extend(m.aveerr[grp])

#toaidx = [(i,x) for i,x in enumerate(alltoa)]
#toaidx.sort(key=lambda y:y[1])

#idx = [i for i,x in toaidx]
#alltoa= np.array(alltoa)
#allres= np.array(allres)
#allerr= np.array(allerr)

#alltoa= alltoa[idx]
#allres= allres[idx]
#allerr= allerr[idx]

#gap = (49375, 50850)
#smallgap = (53150, 53340)
#smallergap = (54220, 54355)

#longest = float(t.end- t.start)/year

#gapsize = gap[1] - gap[0]
#tol = gapsize/longest/year

#saveobject = {'data':np.vstack([alltoa, allres, allerr]), 'span':longest, 'gap':gap, 'smallgap':smallgap, 'smallergap':smallergap, 'tend':t.end, 'tstart':t.start} 

def overlap(binstart, binend, gap, tol):
    binlength = binend - binstart
    if binend <= gap[0]:return True 
    if binstart >= gap[1]:return True 
    ma = max(binstart, gap[0])
    mi = min(binend, gap[1])
    ovlp = mi - ma
    if ovlp/binlength > tol:return False
    return True

saveobject = np.load('timingnoise.npy').item()
data = saveobject['data']
alltoa= data[0,...]
allres= data[1,...]
allerr= data[2,...]

#errorbar(alltoa, allres, yerr=allerr, fmt='.')
#show()

gap = saveobject['gap']
smallgap = saveobject['smallgap']
smallergap = saveobject['smallergap']
gaps=(gap, smallgap, smallergap)

longest = saveobject['span']
tend = saveobject['tend']
tstart = saveobject['tstart']

gapsize = gap[1] - gap[0]
tol = gapsize/longest/year

def SegmentFit(alltoa, allres, allerr, tau, high=longest, shift=0, gaps=(gap, smallgap, smallergap)):
    """
    fit segments of the timing residuals.
    """
    #taus = []
    sigma_z = []
    sigma = lambda tau, p0:1./2./np.sqrt(5)*tau**2*np.sqrt(np.mean(p0**2))

    #tau = np.linspace(np.log10(low), np.log10(longest), N)
    #mjdbins = 10**tau*year
    #mjdbins = [tau*year]
    bin = tau*year
    #for bin in mjdbins:
    pointer = alltoa.size - shift -1
    binend = float(tend) + 0.001
    binstart = binend - bin
    #taus.append(bin/year)
    p0 = []
    while float(tstart) - binstart < 0.1*tol*bin and all([overlap(binstart, binend, g, tol) for g in gaps]):
        segtoa = []
        segres = []
        segerr = []
        while alltoa[pointer] > binstart and pointer >= 0:
            pointer -= 1
            segtoa.append(alltoa[pointer])
            segres.append(allres[pointer])
            segerr.append(allerr[pointer])
        segtoa = np.array(segtoa)
        segres = np.array(segres)
        segerr = np.array(segerr)
        if segtoa.size == 0:
            #print 'no point***'
            pass
        elif segtoa.size <=5:
            pass 
        else:
            x = segtoa - (segtoa.max() - segtoa.min())/2
            x /= year
            y = segres
            w = 1/segerr**2
            pol = np.polyfit(x,y,3, w=w)
            p0.append(pol[0])
            #print 'some points:', segtoa.size, pol[0]
        binend = binstart
        binstart = binend - bin
    return sigma(bin/year, np.array(p0))/1.e6/secperday/year
        #sigma_z.append(sigma(bin/year, np.array(p0)))

    #sigma_z = np.array(sigma_z)/1.e6/secperday/year
    #return taus, sigma_z


sigma_z_10 = SegmentFit(alltoa, allres, allerr, 10)
print sigma_z_10
