import matplotlib 
#matplotlib.use('GTK3cairo')
from datatools.tempo import *
from pylab import *
import matplotlib.pyplot as plt
import numpy as np
from math import *
year = 365.24218967

m = model('1713.Sep.KOM.par')
#m = np.load('1713.Sep.sav')
#m = model('1713.Apr.TN.par')
#t = TOAfile('1713.Nov.Paul.tim')
#m.tempofit(t)
#m.average(groups='allinone')

#wrms = m.wrms['all']
wrms = 0.2124
#wrms = 0.306400344221

CPMSP_sc10= lambda nv, nvdot, T:exp(1.6) * nv**-1.4*abs(nvdot/1.e-15)**1.1*T**2

#Hobbs = lambda nv, nvdot, T:10.**-11.5 * nv*-0.4 * abs(nvdot/1.e-15)**0.8
Hobbs = lambda nv, nvdot, T:10**(-1.37*log10(nv**0.29*abs(nvdot)**-0.55) + 0.52)

Tspan = float(m.FINISH - m.START)/year
nv = float(m.F0[0])
nvdot = float(m.F1[0])

print wrms
print CPMSP_sc10(nv, nvdot, Tspan)
print Hobbs(nv, nvdot, Tspan)

def wrms(res, err):
    nres = [] 
    nerr = [] 
    for i in range(len(err)):
        if not err[i] == 0.:
            nres.append(res[i])
            nerr.append(err[i])
    res = np.array(nres)
    err = np.array(nerr)
    weight = 1./err**2
    wmres = sum(res*weight)/sum(weight)
    mres = np.mean(res)
    sumres = sum(res)
    wsum = sum(weight)
    try: 
        wrms = sqrt(sum(res**2*weight - wmres*sumres)/wsum)
    except ValueError:
        return None 
    return wrms 


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

def SegmentFit(alltoa, allres, allerr, tau, high=longest, shift=0, gaps=(gap, smallgap, smallergap)):
    """
    fit segments of the timing residuals.
    """
    #taus = []
    tol = gapsize/longest/year
    sigma_z = []
    #sigma = lambda tau, p0:1./2./np.sqrt(5)*tau**2*np.sqrt(np.mean(p0**2))

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
    rms = []
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
            pol = np.polyfit(x,y,1, w=w)
            #print pol
            #p0.append(pol[0])
            func = np.poly1d(pol)
            e = y - func(x)
            #rms.append(np.sqrt(np.mean(e**2)))
            rms.append(wrms(e,segerr))

            #print 'some points:', segtoa.size, pol[0]
        binend = binstart
        binstart = binend - bin
    #return sigma(bin/year, np.array(p0))/1.e6/secperday/year
        #sigma_z.append(sigma(bin/year, np.array(p0)))

    #sigma_z = np.array(sigma_z)/1.e6/secperday/year
    #return taus, sigma_z
    #return np.mean(rms), np.std(rms)
    return rms


taus = 10**(np.linspace(np.log10(0.16), np.log10(21), 20))
#print taus
sigma_z_10 = SegmentFit(alltoa, allres, allerr, 10)

sigma_z_list = [SegmentFit(alltoa, allres, allerr, tau) for tau in taus]
t = np.linspace(0.1,30, 30)
RW0 = lambda t:0.1/10.**0.5*t**0.5 + 0.09
RW1 = lambda t:0.1/10.**1.5*t**1.5 + 0.09 
RW2 = lambda t:0.1/10.**2.5*t**2.5 + 0.09
#sigma_z = [s[0] for s in sigma_z_list]
#sigma_err = [s[1] for s in sigma_z_list]
#errorbar(taus, sigma_z, yerr=sigma_err, fmt='.')
ax = plt.subplot(111)
for i, rms in enumerate(sigma_z_list):
    ax.plot([taus[i]]*len(rms), rms, 'k.')
#ax.plot(taus, [0.09]*len(taus), 'k-')
wn, = ax.plot(t, 0.09+t*0, color='k')
rw0, = ax.plot(t, RW0(t), 'k:')
rw1, = ax.plot(t, RW1(t), 'k--')
rw2, = ax.plot(t, RW2(t), 'k-.')
legend([rw2, rw1, rw0, wn], ['RW$_2$', 'RW$_1$', 'RW$_0$', 'WN'], loc=2, numpoints=1)
ax.set_xlim((0.1, 30))
#ax.semilogx()
ax.set_xscale("log")
ax.set_yscale("log")
ax.get_xaxis().set_ticks([0.1, 1, 10])
ax.get_yaxis().set_ticks([0.1, 1, 10])
from matplotlib.ticker import ScalarFormatter, FormatStrFormatter
#majorformator = ScalarFormatter()
majorformator = FormatStrFormatter('%g')
ax.xaxis.set_major_formatter(majorformator)
ax.yaxis.set_major_formatter(majorformator)
ax.set_xlabel('T (yr)', fontsize=14)
ax.set_ylabel('$\sigma^2_{\mathscr{R}, 2}$ ($\mu$s)', fontsize=14)
show()
