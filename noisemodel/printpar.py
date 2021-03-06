from datatools.tempo import *
from TexTable import deluxetable
from round import TexStyle as SF
import numpy as np
from ordereddict import OrderedDict
import decimal
from tools.Coordinate import RA, Dec


secperyear = 3600*24*Decimal('365.24218967')
secperday = 3600 * 24
PI = np.pi

def M1(pf):
    G = float(6.673e-11)
    Msun = float(1.98892e30)
    Tsun = float(4.925490947e-6)
    c = float(2.99792458e8)
    m2 = float(pf.M2[0])
    #I = float(pf.KIN[0])/180*np.pi
    sini = float(pf.SINI[0])
    Pb = float(pf.PB[0])*secperday
    a = float(pf.A1[0])
    #result = sqrt(930.998*m2**3*Pb**2/a**3) - m2
    #return (Pb/2/PI*(sqrt(Tsun*(m2*sin(I))**3/a**3))-m2)
    return (Pb/2/PI*(sqrt(Tsun*(m2*sini)**3/a**3))-m2)

def B(pf):
    return 3.2e19*np.sqrt(np.abs(float(pf.F1[0]/pf.F0[0]**3)))


def Age(pf):
    return np.abs(float(pf.F0[0]/pf.F1[0]/secperyear/2))

def aveDM(pf):
    DMX, DMXErr, DMXR1, DMXR2 = pf.dmxlist
    dmx = [DMX[x] for x in DMX if DMXR2[x]>53200]
    dmxerr = np.std([float(x) for x in dmx])
    return sum(dmx)/len(dmx) , Decimal(str(dmxerr))

def parseerror(val, err):
    s, d ,e = err.as_tuple()
    if int(d[0]) == 1:# and int(d[1]) >= 5:
    #if int(d[0]) == 1 and not int(d[1]) == 0:
    #if int(d[0]) == 1: #and int(d[1]) >= 5:
        errdig = Decimal((s,d[:2],e+len(d)-2))
        #errstr = Decimal((0, d[:2], 0))
        val = str(val.quantize(errdig))
        errstr = ''.join([str(x) for x in err.quantize(errdig).as_tuple()[1]])
    else:
        errdig = Decimal((s,d[:1],e+len(d)-1))
        #errstr = Decimal((0, d[:1], 0))
        val = str(val.quantize(errdig))
        errstr = ''.join([str(x) for x in err.quantize(errdig).as_tuple()[1]])
    if not val.find('E') == -1:
        v,e = val.split('E')
        result = r'$%s(%s)\times10^{%s}$' % (v,errstr, e)
    else:
        result = r'%s(%s)' % (val,errstr)
    return result



Observable = OrderedDict([
        ('RAJ',r'Right Ascension, $\alpha$ (J2000)'),
        ('DECJ',r'Declination, $\delta$ (J2000)'),
        ('F0',r'Spin Frequecy $\nu$~(s$^{-1}$)'),
        ('F1',r"Spin down rate $\nu'$ (s$^{-2}$)"),
        #('F2',r"Second spin frequency derivative $\nu''$ (s$^{-3}$)"),
        ('PMRA',r'Proper motion in $\alpha$, $\nu_{\alpha}=\dot{\alpha}\cos \delta$ (mas~yr$-1$)'),
        ('PMDEC',r'Proper motion in $\delta$, $\nu_{\delta}=\dot{\delta}$ (mas~yr$-1$)'),
        ('PX',r'Parallax, $\pi$ (mas)'),
        ('DM',r'Dispersion Measure# (pc~cm$^{-3}$)'),
        ('SINI', r'$\sin i$, where $i$ is the orbital inclination angle'),
        ('PB',r'Orbital Period, $P_{\rm b}$ (day)'),
        ('E',r'Eccentricity, $e$'),
        ('OM',r'Angle of periastron#, $\omega$ (deg)'),
        ('T0',r'Time of periastron passage, $T_0$ (MJD)'),
        ('PBDOT',r'Change rate of $P_{\rm b}$, $\dot{P}_{\rm b}$ ($10^{-12}$s~s$^{-1}$)'),
        ('A1',r'Projected semi-major axis, $x$ (lt-s)'),
        ('XDOT',r'Change rate of $x$ due to proper motion, $\dot{x}$ (lt-s~s$^{-1}$)'),
        ('M2',r'Companion Mass, $M_c$ ($M_{\odot}$)'),
        ('FD1',r'Profile frequency dependency parameter, FD1 '),
        ('FD2',r'Profile frequency dependency parameter, FD2 '),
        ('FD3',r'Profile frequency dependency parameter, FD3 '),
        ('FD4',r'Profile frequency dependency parameter, FD4 ')
        ])

Fixed = OrderedDict([
    ('EPHEM', r'Solar system ephemeris'),
    ('PEPOCH',r'Refrence epoch for $\alpha$, $\delta$, and $\nu$ (MJD)'),
    ('OMDOT',r'Rate of periastron advance, $\dot{\omega}$ (deg/yr)'),
    ])

Derived = OrderedDict([
        #('KIN',r'Oribital inclination, $i$ (deg)'),
        ('PAASCNODE',r'Position angle of ascending node, $\Omega$ (deg)'),
        ('M1',r'Pulsar mass, $M_{\rm PSR}$ ($M_{\odot}$)'),
        ('B',r'Dipole magnetic field, $B$ (G)'),
        ('Age',r'Characteristic age, $\tau_c$ (yr)'),
        ])


Caption = 'Timing model parameters# from {\it tempo}.'
colnames = ['Parameter', 'EFAC \& EQUAD', 'with jitter model', 'jitter \& red noise model']
parameter = []
data = []

parameter.append(r'\textit{Measured Parameters}')
for k in Observable:
    parameter.append(Observable[k])
parameter.append(r'\textit{Fixed Parameters}')
for k in Fixed:
    parameter.append(Fixed[k])
parameter.append(r'\textit{Derived Parameters}')
for k in Derived:
    parameter.append(Derived[k])

def addmodeldata(m):
    value = []
    value.append('')
    for k in Observable:
        #parameter.append(Observable[k])
        if type(m.__dict__[k]) in (tuple, list):
            val, err = m.__dict__[k]
            if k == 'RAJ':
                ra = RA(m.__dict__[k][0])
                raerr = m.__dict__[k][1]
                value.append('%s:%s:%s' % (ra.HH,ra.MM, parseerror(Decimal(ra.SS),Decimal(raerr))))
            elif k == 'DECJ':
                dec = Dec(m.__dict__[k][0])
                decerr = m.__dict__[k][1]
                value.append('%s:%s:%s' % (dec.dd,dec.mm, parseerror(Decimal(dec.ss),Decimal(decerr))))
            else:
                value.append(parseerror(*m.__dict__[k]))
        else:
            if k == 'DM':
                value.append(parseerror(*aveDM(m)))
            else:
                value.append((m.__dict__[k]))

    value.append('')

    for k in Fixed:
        #parameter.append(Fixed[k])
        try:
            value.append(parseerror(*m.__dict__[k]))
        except:
            if type(m.__dict__[k]) == type(''):
                value.append(m.__dict__[k])
            elif k == 'PEPOCH':
                value.append(m.__dict__[k].quantize(50000))
            elif k == 'OMDOT':
                value.append(m.__dict__[k].quantize(Decimal(0.0001)))
            else:
                value.append(SF(globals()[k](m)))
    #data = [parameter, value]

    value.append('')

    for k in Derived:
        #parameter.append(Derived[k])
        if k == 'PAASCNODE':
            value.append(m.__dict__[k])
        else:
            try:
                value.append(parseerror(*m.__dict__[k]))
            except:
                value.append(SF(globals()[k](m)))
    return value

m = model('normaltempo.par')
value1 = addmodeldata(m)
del m
m = model('glstempo.par')
value2 = addmodeldata(m)
del m
m = model('RNtempo.par')
value3 = addmodeldata(m)
del m

data = [parameter, value1, value2, value3]

comments=[
    #'Numbers in parentheses indicate the uncertainties on the last disgit(s).  Uncertainties on parameters are estimated from a combination of 23 indepedent MCMC chains randomly generated using the Metropolis algorithm. Each chain run for 100000 points after 10000 burn-in steps.', 
    'Numbers in parentheses indicate the uncertainties on the last disgit(s).  Uncertainties on parameters are estimated by the {\it tempo} program using information in the covariance matrix.', 
'   The averaged DM value; See Section 3.2 and Figure 2 for the more discussion.', 
'   See Figure 2 of \citealt{sns+05} for definition.']

#print len(parameter), len(value)

table = deluxetable(Caption=Caption, colsetting='lccc', colnames = colnames, data=data, label="tab:par", comments=comments)

print table

