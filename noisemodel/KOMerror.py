from multiprocessing import Pool, Process, Manager, JoinableQueue, Value
import numpy as np
from datatools.tempo import *
from decimal import *
from tempfile import mkdtemp
from threadit import threadit
import os, sys
from pylab import *

parfile = '1713.Sep.KOM.par'
toafile = '1713.Sep.T2.tim'

def calchisq(offset):
    cwd=os.getcwd()
    tmpdir = cwd+'/.'+uniquename()
    if not tmpdir == None:
        if os.path.exists(tmpdir):
            os.chdir(tmpdir)
        else:
            os.mkdir(tmpdir)
            os.chdir(tmpdir)
    os.system('cp %s/%s %s/%s' % (cwd, parfile, tmpdir, parfile))
    m = model(parfile)
    m.PAASCNODE += Decimal(offset)
    m.freezeall()
    m.write(parfile)
    t = TOAfile('../' + toafile)
    m.tempofit(t, GLS=True)
    os.chdir(cwd)
    return m.chisq

#print calchisq(0.)

offsets = [ [x] for x in np.linspace(-10,10,20)]
chisqs = threadit(calchisq, offsets)
np.save('KOMs', (offsets,chisqs))
plot(offsets, chisqs, 'o')
show()
