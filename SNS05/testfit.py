from pylab import *
from datatools.tempo import *

md = model('1713.sns.par')
tf = TOAfile('1713.sns.tim')
#md.tempofit(tf, GLS=True)
md.tempo2fit(tf)
md.average()

#a1 = subplot(211)
#a2 = subplot(212)


#md.plot('date', 'DMX', ax=a1)
#md.plot('date', 'averes', ax=a2, LegendOn=True)
md.plot('date', 'averes',  LegendOn=True)

show()
