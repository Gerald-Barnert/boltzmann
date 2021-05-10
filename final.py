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
        pars.set_matter_power(redshifts=zz, kmax=2.0)
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

#pars, results = camb_code[0], camb_code[1]
#results = camb_code[1]
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

def growth_rate(z, run, code):
    if code == "camb":
        pars = run[0]
        results = run[1]
        #results = camb.get_results(pars)
        #pars.set_matter_power(redshifts=z, kmax=2.0)
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

def Angular_power_spectrum(z, run, code, l_max, units):
    if code == "camb":
        pars = run[0]
        results = run[1]
        pars.InitPower.set_params()
        pars.set_for_lmax(l_max, lens_potential_accuracy=0)
        powers =results.get_cmb_power_spectra(pars, CMB_unit= units)
        totCl=powers['total']
        ls = np.arange(totCl.shape[0])
        Cl_TT = totCl[:,0]
        #Cl_EE = totCl[:,1]
        #Cl_BB = totCl[:,2]
        #Cl_ET = totCl[:,3]
        return ls[3:], Cl_TT[3:]

    if code == "class":
        Cl = run.lensed_cl(l_max)
        ls = Cl['ell']
        if units == "muK":
            Cl_TT = Cl['tt']*ls*(ls+1) / (2*np.pi) * (run.T_cmb()*10**6)**2
        if units == None:
            Cl_TT = Cl['tt']*ls*(ls+1) / (2*np.pi)
        #Cl_EE = Cl['ee']*ls*(ls+1) / (2*np.pi)
        #Cl_BB = Cl['bb']*ls*(ls+1) / (2*np.pi)
        #Cl_ET = Cl['te']*ls*(ls+1) / (2*np.pi)
        return ls[3:], Cl_TT[3:]

def Matter_power_spectrum(z, run, code, l_max, kmin=1e-4, kmax=2, number_points=2000):
    if code == "camb":
        pars = run[0]
        results = run[1]
        pars.InitPower.set_params()
        pars.set_matter_power(redshifts=z, kmax=kmax)
        pars.set_for_lmax(l_max, lens_potential_accuracy=0)
        kh,z, pk = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = number_points)
        P = [] 
        for i in range(0, len(z)):
            P.append(pk[i,:])
        return kh, P

    if code == "class":
        kk = np.linspace(kmin, kmax, number_points)
        Pk_z = []
        P = []
        h = run.h()
        print(h)
        for zz in z:
            Pk = []
            for k in kk:
                Pk.append(run.pk(k*h,zz) * h**3)
            Pk_z.append(Pk)
        #for i in range(0, len(z)):
        #    P.append(Pk_z[i])
        return kk, Pk_z

p1 = Matter_power_spectrum(zz, camb_code, "camb", l_max=2500)
p2 = Matter_power_spectrum(zz, class_code, "class", l_max=2500)

a1 = Angular_power_spectrum(zz, camb_code, "camb", l_max=2500, units=None)
a2 = Angular_power_spectrum(zz, class_code, "class", l_max=2500, units=None)

huble = HubbleParameter(zz, camb_code, "camb")
huble2 = HubbleParameter(zz, class_code, "class")

d_a = Angular_diameter_distance(zz, camb_code, "camb")
d_a2 = Angular_diameter_distance(zz, class_code, "class")

f1 = growth_rate(zz, camb_code, "camb")
f2 = growth_rate(zz, class_code, "class")

#plt.loglog(p1[0], p1[1][0])
#plt.loglog(p2[0], p2[1][0])

#plt.plot(a1[0], a1[1])
#plt.plot(a2[0], a2[1])

#plt.plot(zz,d_a)
#plt.plot(zz,d_a2)

#plt.plot(zz, huble)
#plt.plot(zz, huble2)
plt.show()

