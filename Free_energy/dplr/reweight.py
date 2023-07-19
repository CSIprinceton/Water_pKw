import numpy as np
import matplotlib.pyplot as plt
from sys import argv

colvar=np.loadtxt('COLVAR').transpose()

Kb=0.00831446 #kJ/mol
T=330 #K
L=float(argv[1])
rho=0.0334

nbins=int(np.max(colvar[1])/0.05)
bias=colvar[3]+colvar[5]
hist,bins=np.histogram(colvar[1],bins=nbins,weights=np.exp(bias/(Kb*T)))
bins=(bins[1:]+bins[:-1])/2.
hist/=4*np.pi*bins**2*(bins[1]-bins[0])
hist/=4*np.pi*np.trapz(bins**2*hist,bins)*(bins[1]-bins[0])/rho

fes=-np.log(hist)

filt=bins<L/2.
np.savetxt('fes.dat',np.vstack((bins[filt],fes[filt])).transpose())

