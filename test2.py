# -*- coding: utf-8 -*-
import sys
reload(sys)
sys.setdefaultencoding('utf8')

import camb
import classy
from classy import Class
from camb import model
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

zz = np.arange(0,3,0.05)
pars = camb.CAMBparams()
pars.set_cosmology(H0=70, ombh2=0.022, omch2=0.122)
#pars.InitPower.set_params(ns=0.965)
pars.set_matter_power(redshifts=zz, kmax=2.0)
results = camb.get_results(pars)
kh,z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints = 2000)

k = kh / pars.h
print(k)



omb = pars.ombh2
omc = pars.omch2

print((omb + omc)/(pars.h * pars.h))
print(pars.h)
print(pars.H0)

H = []
D_A = []
F_AP = []
A = results.get_BAO(zz, pars)
for i in range(0, len(zz)):
    hubble = A[i][1]
    DA = A[i][2]
    F = A[i][3]
    H.append(hubble)
    D_A.append(DA)
    F_AP.append(F)

fsigma8 = results.get_fsigma8()
sigma8 = results.get_sigma8()
f_growth_rate = np.flip(fsigma8/sigma8)
print(sigma8)

plt.plot(zz, f_growth_rate)
plt.ylabel(r'Growth rate $f_g(z)$')
plt.xlabel(r'$z$')
plt.savefig('plots/growth_rate_f(z)_camb')
plt.show()

pars.InitPower.set_params(As=2e-9, ns=0.965, r=0)
pars.set_for_lmax(2500, lens_potential_accuracy=0)
powers =results.get_cmb_power_spectra(pars, CMB_unit= None)
for name in powers: print(name)
totCl=powers['total']
ls = np.arange(totCl.shape[0])
Cl_TT = totCl[:,0]
Cl_EE = totCl[:,1]
Cl_BB = totCl[:,2]
Cl_ET = totCl[:,3]

'''
plt.plot(ls[3:], Cl_TT[3:])
plt.ylabel(r'$l(l+1)\:C^{TT}_l\: / 2\pi$')
plt.xlabel(r'$Multipole\: moment\: l$')
plt.show()
'''
#plt.plot(ls[3:], Cl_EE[3:])
#plt.show()
#plt.plot(ls[3:], Cl_BB[3:])
#plt.show()
#plt.plot(ls[3:], Cl_ET[3:])
#plt.show()


'''
for i in range(0, len(zz), 10):
    plt.loglog(kh, pk[i,:], label = '{}'.format(zz[i]))

plt.legend(title='z')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$')
plt.xlabel(r'$k\:[h/Mpc]$')
#plt.savefig('plots/P(k)_camb')
plt.show()

plt.plot(zz,H)
plt.ylabel(r'$H(z)\:[Km\:s⁻¹\:Mpc⁻¹]$')
plt.xlabel('z')
#plt.savefig('plots/H(z)_camb')
plt.show()

plt.plot(zz, D_A)
plt.ylabel(r'$d_A(z)\: [Mpc]$')
plt.xlabel('z')
#plt.savefig('plots/d_A_camb')
plt.show()
'''


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