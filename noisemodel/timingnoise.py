from datatools.tempo import *
from pylab import *
import numpy as np
import scipy.optimize as opt

year = 365.24218967

#t = TOAfile('1713.Oct.t2.toa')
#m = model('1713.Oct.jump.par')
#t = TOAfile('1713.Nov.Paul.tim')
#t = TOAfile('1713.sns.tim')
#m = model('1713.Nov.test.par')
#m = model('1713.Nov.test.par')
#m = model('1713.Jan.mcmc.par')
#m = model('1713.Mar.mcmc.par')
#m = model('1713.Mar.TN.par')
#m = model('1713.dmx.mcmc.par')
#m = model('1713.sns.par')
#m.tempofit(t)
#mx = model(m.newpar.parfile)
#ratio = mx.F2[0]/mx.F0[0]
#TZRMJD = mx.TZRMJD
#resfunc = lambda x:float(ratio)*float(((Decimal(x) - TZRMJD)*secperday)**3) * -1.e6/6
#mx = model('1713.Nov.DMX.par')

#mx.freezeall()
#mx.F2[0] = Decimal(0)
#mx.F3[0] = Decimal(0)
#mx.F4[0] = Decimal(0)
#mx.F5[0] = Decimal(0)
#mx.F6[0] = Decimal(0)
#mx.F7[0] = Decimal(0)
#mx.F8[0] = Decimal(0)
#mx.thawall('JUMP')
#mx.write('1713.Mar.TN.par')
#del m, mx, t
#mx = model('1713.Mar.TN.par')
mx = model('1713.Sep.KOM.par')
t = TOAfile('1713.Sep.T2.tim')
mx.tempofit(t, GLS=True)
mx.average(lapse = 0.25)
#mx.plot('date', 'averes')
#m.plot('date', 'prefit')
#show()
groups = mx.averes.keys()
groups.sort()
for grp in groups:
    errorbar(mx.avetoa[grp], mx.averes[grp], yerr=mx.aveerr[grp], fmt = '.', mew=0, label=grp)

np.save(open('1713.Sep.sav', 'wb'), mx)

#def TNfunc(DC):
    #res = []
    #for grp in groups:
        #res.extend([(mx.averes[grp][i] - resfunc(x) - DC)**2 for i,x in enumerate(mx.avetoa[grp])])
    #return np.sum(res)

#result = opt.fmin(TNfunc, 0.)
##print result

#mjds = np.arange(mx.START, mx.FINISH, 100)
#TNsig = [resfunc(x) + result[0] for x in mjds]

#plot(mjds, TNsig, '--')
xlabel(r'MJD (day)')
ylabel(r'Residual ($\mu$s)')
show()

alltoa = []
allres = []
allerr = []
for grp in groups:
    alltoa.extend(mx.avetoa[grp])
    allres.extend(mx.averes[grp])
    allerr.extend(mx.aveerr[grp])

toaidx = [(i,x) for i,x in enumerate(alltoa)]
#newtoaidx = sorted(toaidx, key=lambda y:y[1])
toaidx.sort(key=lambda y:y[1])

idx = [i for i,x in toaidx]
alltoa= np.array(alltoa)
allres= np.array(allres)
allerr= np.array(allerr)

alltoa= alltoa[idx]
allres= allres[idx]
allerr= allerr[idx]

#errorbar(alltoa, allres, yerr=allerr, fmt = '.', mew=0)
#show()

gap = (49375, 50850)
smallgap = (53150, 53340)
smallergap = (54220, 54355)

longest = float(t.end- t.start)/year

gapsize = gap[1] - gap[0]
tol = gapsize/longest/year

saveobject = {'data':np.vstack([alltoa, allres, allerr]), 'span':longest, 'gap':gap, 'smallgap':smallgap, 'smallergap':smallergap, 'tend':t.end, 'tstart':t.start} 


def overlap(binstart, binend, gap, tol):
    binlength = binend - binstart
    if binend <= gap[0]:return True 
    if binstart >= gap[1]:return True 
    ma = max(binstart, gap[0])
    mi = min(binend, gap[1])
    ovlp = mi - ma
    if ovlp/binlength > tol:return False
    return True

tau = np.linspace(np.log10(0.25), np.log10(longest), 10)
mjdbins = 10**tau*year
print tol

#saveobject = {}

for bin in mjdbins:
    #saveobject[bin] = []
    #print bin
    pointer = alltoa.size - 1
    binend = float(t.end) + 0.0001
    binstart = binend - bin
    #while float(t.start) - binstart < 0.1*tol*bin and overlap(binstart, binend, gap, tol) and overlap(binstart, binend, smallgap, tol) and overlap(binstart, binend, smallergap, tol):
    while float(t.start) - binstart < 0.1*tol*bin and overlap(binstart, binend, gap, tol)  and overlap(binstart, binend, smallergap, tol):
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
        #saveobject[bin].append([segtoa, segres, segerr])
        #print segtoa, bin, binstart, binend, gap
        x = segtoa - (segtoa.max() - segtoa.min())/2
        y = segres
        w = 1/segerr**2
        #if len(x) == 0: print bin, pointer 
        #pol = np.polyfit(x,y,3)
        #print pol
        #polyfunc = np.poly1d(pol)
        #errorbar(segtoa, segres, segerr, fmt='.', mew=0)
        #plot(segtoa, polyfunc(x), '--')
        binend = binstart
        binstart = binend - bin

np.save('timingnoise', saveobject)
#show()
