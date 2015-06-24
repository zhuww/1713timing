from pylab import *
from datatools.tempo import *


md = model('1713.sns.par')
for par in ['F2', 'F3', 'F4']:
    md.__dict__[par][0] = Decimal(0)

md.freezeall()
md.thawall('F0')
md.thawall('F1')

tf = TOAfile('1713.sns.tim')



md.tempo2fit(tf)
md.plot('mjd', 'res')

show()
