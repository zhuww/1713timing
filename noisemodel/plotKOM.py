import numpy as np
from pylab import *

KOM, chisq = np.load('KOMs.npy')
kom = np.array([x[0] for x in KOM])

p = np.polyfit(kom, chisq, 2)
print  p
poly = lambda x:p[0] * x**2 + p[1] * x + p[2]

plot(kom, chisq, 'o')
plot(kom, poly(kom), '--')
show()

a = p[0]
b = (0.-p[1])/2./a
#dx = np.sqrt(2.3/a)
dx = np.sqrt(1./a)
print b ,dx
