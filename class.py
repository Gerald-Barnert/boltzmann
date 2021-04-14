# -*- coding: utf-8 -*-
import sys
reload(sys)
sys.setdefaultencoding('utf8')

from classy import Class
import numpy as np
import matplotlib.pyplot as plt

LambdaCDM = Class()
LambdaCDM.set({'omega_b': 0.022, 'omega_cdm': 0.122, 'h': 0.70, 'z_max_pk':4})
LambdaCDM.set({'output': 'tCl, lCl, pCl, mPk' , 'lensing':'yes', 'P_k_max_1/Mpc': 3.0 })
LambdaCDM.compute()

kk = np.linspace(1e-4, 2, 2000)
Pk = []
Pk_z = []
h = LambdaCDM.h()
zz = np.arange(0,3.5,0.05)

for z in zz:
    Pk = []
    for k in kk:
        Pk.append(LambdaCDM.pk(k*h,z) * h**3)
    Pk_z.append(Pk)

for i in range(0, len(zz)):
    plt.loglog(kk, Pk_z[i])

plt.ylabel(r'$p(k)\:[(Mpc/s)³]$')
plt.xlabel(r'$k\:[h/Mpc]$')
plt.show()
#plt.savefig('p(k)_class')
plt.clf()

c = 3e5

H = []
D_A = []
for z in zz:
    H.append(LambdaCDM.Hubble(z)*c)
    D_A.append(LambdaCDM.angular_distance(z))


print(H[0])
plt.plot(zz,H)
plt.ylabel(r'$H(z)\:[Km\:s⁻¹\:Mpc⁻¹]$')
plt.xlabel('z')
plt.show()

plt.plot(zz, D_A)
plt.ylabel(r'$d_A\:(z) [Mpc]$')
plt.xlabel('z')
plt.show()
