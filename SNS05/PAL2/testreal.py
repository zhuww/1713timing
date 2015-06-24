from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import PAL2.PALutils as PALutils
import PAL2.PALmodels as PALmodels
from datatools.MJD import MJD_to_datetime
#import PALmodels
#import PALutils
def dailyAve(times, res, err, ecorr=None, dt=10, flags=None):

    isort = np.argsort(times)
    
    bucket_ref = [times[isort[0]]]
    bucket_ind = [[isort[0]]]
    
    for i in isort[1:]:
        if times[i] - bucket_ref[-1] < dt:
            bucket_ind[-1].append(i)
        else:
            bucket_ref.append(times[i])
            bucket_ind.append([i])
    
    avetoas = np.array([np.mean(times[l]) for l in bucket_ind],'d')
    if flags is not None:
        aveflags = np.array([flags[l[0]] for l in bucket_ind])

    aveerr = np.zeros(len(bucket_ind))
    averes = np.zeros(len(bucket_ind))

    for i,l in enumerate(bucket_ind):
        M = np.ones(len(l))
        C = np.diag(err[l]**2) 
        if ecorr is not None:
            C += np.ones((len(l), len(l))) * ecorr[l[0]]

        avr = 1/np.dot(M, np.dot(np.linalg.inv(C), M))
        aveerr[i] = np.sqrt(avr)
        averes[i] = avr * np.dot(M, np.dot(np.linalg.inv(C), res[l]))
 
        
    if flags is not None:
        return avetoas, averes, aveerr, aveflags
    else:
        return avetoas, aveerr, averes

h5file = './sns.h5'
pulsar = ['J1713+0747']

model = PALmodels.PTAmodels(h5file, pulsars=pulsar)
fullmodel = model.makeModelDict(incRedNoise=True, noiseModel='powerlaw',
                                separateEfacsByFreq=True, separateEquadsByFreq=True,
                                separateJitterEquadByFreq=True, incEquad=True,
                                incJitterEquad=False, likfunc='mark6')

for sig in fullmodel['signals']:
    if sig['stype'] == 'efac':
        #print sig['flagvalue']
        #if np.any([e in sig['flagvalue'] for e in ['M4', 'M3', 'ABPP']]):
        #if True:
        sig['bvary'] = [False]
model.initModel(fullmodel, memsave=True, fromFile=False,
                verbose=True, write='no')
chain = np.loadtxt('chain_1.txt')
#burn = 0.25 * chain.shape[0]
maxpars = chain[np.argmax(chain[:,-3]), :-4]
#maxpars = np.array([ -8.35925582, -9.32412922, -6.89705264,
                    #-6.78209842,-6.8684193,-8.50634118,
                    #-14.85146377,5.43897386])
res_red, err_red = model.create_realization(maxpars, signal='red', incJitter=False)
res_tm, err_tm = model.create_realization(maxpars, signal='tm', incJitter=False)
whiteres = model.psr[0].residuals - res_tm[0] - res_red[0]

#ecorr = np.dot(model.psr[0].Ttmat[:,-njitter:], model.psr[0].Qamp)
#avetoas, aveerr, averes = dailyAve(model.psr[0].toas, whiteres, np.sqrt(model.psr[0].Nvec), ecorr=ecorr)

#plt.errorbar(model.psr[0].toas, model.psr[0].residuals*1e6, model.psr[0].toaerrs*1e6, fmt='.')
#plt.errorbar(model.psr[0].toas, whiteres*1e6, np.sqrt(model.psr[0].Nvec)*1e6, fmt='.')
ind = np.argsort(model.psr[0].toas)
plt.plot(model.psr[0].toas[ind], res_red[0][ind]*1e6, lw=2, ls='--', color='k')
plt.fill_between(model.psr[0].toas[ind], (res_red[0]-err_red[0])[ind]*1e6,
                         (res_red[0]+err_red[0])[ind]*1e6, color='gray', alpha=0.5)

#red, = plt.plot(model.psr[0].toas, res_red[0]*1e6, '.')
#tm, = plt.plot(model.psr[0].toas, res_tm[0]*1e6, '.')

#plt.errorbar(model.psr[0].toas, whiteres*1.e6, np.sqrt(model.psr[0].Nvec)*1e6, fmt='.')
plt.errorbar( model.psr[0].toas, (model.psr[0].residuals - res_tm[0])*1e6 , np.sqrt(model.psr[0].Nvec)*1e6, fmt='.')
plt.xlabel('MJD')
plt.ylabel('residual ($\mu$s)')
#plt.legend([red, tm], ['red', 'tm'])
plt.show()
