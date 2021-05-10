# -*- coding: utf-8 -*-
import camb
import classy
from classy import Class
from camb import model
import numpy as np
import matplotlib.pyplot as plt

c = 3e5

def run(code, Hubble0, om_b, om_cdm):
    if code == "camb":
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=Hubble0, ombh2=om_b, omch2=om_cdm)
        results = camb.get_results(pars)
        return pars, results
    
    if code == "class":
        LambdaCDM = Class()
        h = 0.7
        LambdaCDM.set({'omega_b': 0.022, 'omega_cdm': 0.122, 'h': h, 'z_max_pk':4})
        LambdaCDM.set({'output': 'tCl, lCl, pCl, mPk','lensing':'yes', 'P_k_max_1/Mpc': 3.0 })
        LambdaCDM.compute()
        return LambdaCDM

zz = np.arange(0,3,0.05)
camb_code = run("camb", 70, 0.022, 0.122)
class_code = run("class", 70, 0.022, 0.122)

pars, results = camb_code[0], camb_code[1]
results = camb_code[1]
#print(results)
#pars = camb.CAMBparams()
#pars.set_cosmology(H0=70, ombh2=0.022, omch2=0.122)
#pars.InitPower.set_params(ns=0.965)
#pars.set_matter_power(redshifts=zz, kmax=2.0)
#results = camb.get_results(pars)
#print(results)
#kh,z, pk = results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints = 2000)
#fsigma8 = results.get_fsigma8()


def HubbleParameter(z, run, code):
    if code == "camb":
        pars, results = run[0], run[1]
        H = []
        A = results.get_BAO(zz, pars)
        for i in range(0, len(zz)):
            hubble = A[i][1]
            H.append(hubble)
        return H

    if code == "class":
        H = []
        for z in zz:
            H.append(run.Hubble(z)*c)
        return H

def Angular_diameter_distance(z, run, code):
    if code == "camb":
        pars, results = run[0], run[1]
        D_A = []
        A = results.get_BAO(zz, pars)
        for i in range(0, len(zz)):
            DA = A[i][2]
            D_A.append(DA)
        return D_A

    if code == "class":
        D_A = []
        for z in zz:
            D_A.append(run.angular_distance(z))
        return D_A

def growth_rate(z, code):
    if code == "camb":
        run1 = run("camb", 70, 0.022, 0.122)
        pars = run1[0]
        results = run1[1]
        results = camb.get_results(pars)
        pars.set_matter_power(redshifts=z, kmax=2.0)
        results.get_matter_power_spectrum(minkh=1e-4, maxkh=2, npoints = 2000)
        fsigma8 = results.get_fsigma8()
        sigma8 = results.get_sigma8()
        f_growth_rate = np.flip(fsigma8/sigma8)
        return f_growth_rate

    if code == "class":
        growth_rate_f = []
        for z in zz:
            growth_rate_f.append(run.scale_independent_growth_factor_f(z))
        return growth_rate_f

huble = HubbleParameter(zz, camb_code, "camb")
huble2 = HubbleParameter(zz, class_code, "class")

d_a = Angular_diameter_distance(zz, camb_code, "camb")
d_a2 = Angular_diameter_distance(zz, class_code, "class")

f1 = growth_rate(zz, "camb")
f2 = growth_rate(zz, "class")

plt.plot(zz, f1)
plt.plot(zz, f2)


#plt.plot(zz,d_a)
#plt.plot(zz,d_a2)

#plt.plot(zz, huble)
#plt.plot(zz, huble2)
plt.show()

