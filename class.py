# -*- coding: utf-8 -*-
import sys
reload(sys)
sys.setdefaultencoding('utf8')
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

LambdaCDM = Class()
LambdaCDM.set({'omega_b': 0.022, 'omega_cdm': 0.122, 'h': 0.70, 'z_max_pk':4})
LambdaCDM.set({'output': 'tCl, lCl, pCl, mPk','lensing':'yes', 'P_k_max_1/Mpc': 3.0 })
LambdaCDM.compute()

print(LambdaCDM.Omega0_cdm())
print(LambdaCDM.h())
print(0.122 / (LambdaCDM.h() * LambdaCDM.h()))

kk = np.linspace(1e-4, 2, 2000)
Pk = []
Pk_z = []
h = LambdaCDM.h()
zz = np.arange(0,3,0.05)

for z in zz:
    Pk = []
    for k in kk:
        Pk.append(LambdaCDM.pk(k*h,z) * h**3)
    Pk_z.append(Pk)

c = 3e5

H = []
D_A = []
growth_rate_f = []
for z in zz:
    H.append(LambdaCDM.Hubble(z)*c)
    D_A.append(LambdaCDM.angular_distance(z))
    growth_rate_f.append(LambdaCDM.scale_independent_growth_factor_f(z))

'''
plt.plot(zz, growth_rate_f)
plt.ylabel(r'Growth rate $f_g(z)$')
plt.xlabel(r'$z$')
#plt.savefig('plots/growth_rate_f(z)_class')
plt.show()
'''

Cl = LambdaCDM.lensed_cl(2500)
print(Cl.viewkeys())

ls = Cl['ell']
Cl_TT = Cl['tt']*ls*(ls+1) / (2*np.pi) * (LambdaCDM.T_cmb()*10**6)**2
Cl_EE = Cl['ee']*ls*(ls+1) / (2*np.pi)
Cl_BB = Cl['bb']*ls*(ls+1) / (2*np.pi)
Cl_ET = Cl['te']*ls*(ls+1) / (2*np.pi)


plt.plot(ls[3:], Cl_TT[3:])
plt.ylabel(r'$l(l+1)\:C^{TT}_l\: / 2\pi$')
plt.xlabel('Multipole moment l')
#plt.savefig('plots/C(l)_class')
plt.show()
#plt.plot(ls[3:], Cl_EE[3:])
#plt.show()
#plt.plot(ls[3:], Cl_BB[3:])
#plt.show()
#plt.plot(ls[3:], Cl_ET[3:])
#plt.show()



'''
for i in range(0, len(zz), 10):
    plt.loglog(kk, Pk_z[i], label='{}'.format(zz[i]))

plt.legend(title='z')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$')
plt.xlabel(r'$k\:[h/Mpc]$')
#plt.savefig('plots/P(k)_class')
plt.show()
#plt.savefig('p(k)_class')
plt.clf()

print(H[0])
plt.plot(zz,H)
plt.ylabel(r'$H(z)\:[Km\:s⁻¹\:Mpc⁻¹]$')
plt.xlabel('z')
#plt.savefig('plots/H(z)_class')
plt.show()

plt.plot(zz, D_A)
plt.ylabel(r'$d_A(z)\: [Mpc]$')
plt.xlabel('z')
#plt.savefig('plots/d_A_class')
plt.show()
'''