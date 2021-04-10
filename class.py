from classy import Class
import numpy as np
import matplotlib.pyplot as plt

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