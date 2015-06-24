import numpy as np

data = np.genfromtxt('chain_1.txt', dtype='f')
pars = np.genfromtxt('pars.txt', dtype='|S30')

#print pars.size, pars.shape
#print data.size, data.shape

j = data[:,-3].argmax()
print j
for i in range(pars.size):
    par = str(pars[i])
    if par.startswith('efac'):
        print pars[i].replace('efac-', 'T2EFAC -i '), data[j,i]
    elif par.startswith('equad'):
        print pars[i].replace('equad-', 'T2EQUAD -i '), 10.**(data[j,i] + 6)
    elif par.startswith('RN-Amplitude'):
        print 'TNRedAmp', data[j,i]
    elif par.startswith('RN-spectral-index'):
        print 'TNRedGam', data[j,i]
