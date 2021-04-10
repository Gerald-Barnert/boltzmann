import camb
import classy
from classy import Class
from camb import model
import numpy as np
import matplotlib.pyplot as plt

'''
pars = camb.CAMBparams()
pars.set_cosmology(H0=70, ombh2=0.022, omch2=0.122)
#pars.InitPower.set_params(ns=0.965)
pars.set_matter_power(redshifts=[1,2,3], kmax=1.0)
results = camb.get_results(pars)
kh,z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints = 200)

plt.loglog(kh, pk[0,:])
plt.loglog(kh, pk[1,:])
plt.loglog(kh, pk[2,:])
plt.show()
'''

LambdaCDM = Class()
LambdaCDM.set({'omega_b': 0.022, 'omega_cdm': 0.12, 'h': 0.73, 'z_max_pk':4})
LambdaCDM.set({'output': 'tCl, lCl, pCl, mPk' , 'lensing':'yes', 'P_k_max_1/Mpc': 3.0 })
LambdaCDM.compute()

kk = np.linspace(1e-4, 2, 2000)
Pk = []
Pk_z = []
h = LambdaCDM.h()
zz = np.arange(0,3,1)
print(zz)

for z in zz:
    Pk = []
    for k in kk:
        Pk.append(LambdaCDM.pk(k*h,z) * h**3)
    Pk_z.append(Pk)

plt.loglog(kk, Pk_z[0])
plt.loglog(kk, Pk_z[1])
plt.loglog(kk, Pk_z[2])

plt.show()
