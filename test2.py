import camb
import classy
from classy import Class
from camb import model
import numpy as np
import matplotlib.pyplot as plt

pars = camb.CAMBparams()
pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
#pars.InitPower.set_params(ns=0.965)
pars.set_matter_power(redshifts=[1,2,3], kmax=1.0)
results = camb.get_results(pars)
kh,z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints = 200)

plt.loglog(kh, pk[0,:])
plt.loglog(kh, pk[1,:])
plt.loglog(kh, pk[2,:])
plt.show()

