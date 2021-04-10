import camb
import classy
from classy import Class
from camb import model
import numpy as np
import matplotlib.pyplot as plt

zz = np.arange(0,3.5,0.05)
pars = camb.CAMBparams()
pars.set_cosmology(H0=70, ombh2=0.022, omch2=0.122)
#pars.InitPower.set_params(ns=0.965)
pars.set_matter_power(redshifts=zz, kmax=2.0)
results = camb.get_results(pars)
kh,z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints = 2000)

for i in range(0, len(zz)):
    plt.loglog(kh, pk[i,:])
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
zz = np.arange(0,3.5,0.05)
print(zz)

for z in zz:
    Pk = []
    for k in kk:
        Pk.append(LambdaCDM.pk(k*h,z) * h**3)
    Pk_z.append(Pk)

for i in range(0, len(zz)):
    plt.loglog(kk, Pk_z[i])

plt.show()
'''