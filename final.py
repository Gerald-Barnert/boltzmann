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

c = 3e5

def run(code, Hubble0, Om_b, Om_cdm, As=2.3e-9, ns=0.96,r=0, z_pk_max=4,tau=0.09, k_max=3):
    if code == "camb":
        h = float(Hubble0) /100
        om_b = Om_b*h**2
        om_cdm = Om_cdm*h**2
        pars = camb.CAMBparams()
        pars.set_cosmology(H0=Hubble0, ombh2=om_b, omch2=om_cdm, tau=tau)
        results = camb.get_results(pars)
        pars.InitPower.set_params(As=As, ns=ns, r=r)
        pars.set_matter_power(redshifts=zz, kmax=k_max)
        return pars, results
    
    if code == "class":
        LambdaCDM = Class()
        h = float(Hubble0) /100
        om_b = Om_b*h**2
        om_cdm = Om_cdm*h**2
        LambdaCDM.set({'omega_b': om_b, 'omega_cdm': om_cdm, 'h': h, 'z_max_pk':z_pk_max, 'A_s':As, 'n_s':ns, 'tau_reio':tau})
        LambdaCDM.set({'output': 'tCl, lCl, pCl, mPk','lensing':'yes', 'P_k_max_1/Mpc': k_max })
        LambdaCDM.compute()
        return LambdaCDM

def HubbleParameter(z, run, code):
    if code == "camb":
        pars, results = run[0], run[1]
        A = results.get_BAO(zz, pars)
        H = A[:,1] 
        return H

    if code == "class":
        hubble_vect = np.vectorize(run.Hubble)
        H = hubble_vect(zz)*c
        return H

def Angular_diameter_distance(z, run, code):
    if code == "camb":
        pars, results = run[0], run[1]
        A = results.get_BAO(zz, pars)
        D_A = A[:,2]
        return D_A

    if code == "class":
        DA_vect = np.vectorize(run.angular_distance)
        D_A = DA_vect(zz)
        return D_A

def growth_rate(z, run, code, min_k=1e-4, max_k=3, n_points = 2000):
    if code == "camb":
        pars, results = run[0], run[1]
        results.get_matter_power_spectrum(minkh=min_k, maxkh=max_k, npoints = n_points)
        fsigma8, sigma8 = results.get_fsigma8(), results.get_sigma8()
        f_growth_rate = np.flip(fsigma8/sigma8)
        return f_growth_rate

    if code == "class":
        f_vect = np.vectorize(run.scale_independent_growth_factor_f)
        growth_rate_f = f_vect(zz)
        return growth_rate_f

def Angular_power_spectrum(z, run, code, l_max, units, lens_accuracy=0):
    if code == "camb":
        pars, results = run[0], run[1]
        pars.set_for_lmax(l_max, lens_potential_accuracy=lens_accuracy)
        powers =results.get_cmb_power_spectra(pars, CMB_unit= units)
        totCl=powers['total']
        results.get_cmb_unlensed_scalar_array_dict
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


def P_k_camb_plot(z, pk, zplot, P_z):
    for i in range(0, len(z)):
        for z_j in zplot:
            if z[i] == z_j:
                P_z.append(pk[i,:])

def P_k_class(z, kk, h, zplot, Pk, Pk_z, run, pk_vect, plot):
    for zz in z:
        P = pk_vect(kk*h,zz) * h**3
        Pk.append(P)
        if plot == 'yes':
            for z_j in zplot:
                if zz == z_j:
                    Pk_zz = pk_vect(kk*h,z_j) * h**3
                    Pk_z.append(Pk_zz)


def Matter_power_spectrum(z, run, code, l_max, zplot, kmin=1e-4, kmax=2, number_points=2000, plot="no"):
    if code == "camb":
        pars = run[0]
        results = camb.get_results(pars)
        kh,z, pk = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = number_points)
        P_z = []
        if plot == "yes":
            P_k_camb_plot(z, pk, zplot, P_z)
            return kh, pk, P_z
        if plot == "no":
            return kh, pk

    if code == "class":
        kk = np.linspace(kmin, kmax, number_points)
        Pk = []
        P = []
        Pk_z = []
        h = run.h()
        pk_vect = np.vectorize(run.pk)
        P_k_class(z, kk, h, zplot, Pk, Pk_z, run, pk_vect, plot)
        if plot == 'yes':
            return kk, Pk, Pk_z
        if plot == 'no':
            return kk, Pk

z_plot = [0,1,2]
zmin = 0
zmax = 3
step = 0.05
zz = np.arange(zmin,zmax,step)
#camb_code = run("camb", 70, 0.022, 0.122)
#class_code = run("class", 70, 0.022, 0.122)

camb_code = run("camb", 70, 0.04, 0.24)
class_code = run("class", 70, 0.04, 0.24)

p1 = Matter_power_spectrum(zz, camb_code, "camb", l_max=2500, zplot=z_plot, plot='yes')
p2 = Matter_power_spectrum(zz, class_code, "class", l_max=2500, zplot=z_plot, plot='yes')

a1 = Angular_power_spectrum(zz, camb_code, "camb", l_max=2500, units=None)
a2 = Angular_power_spectrum(zz, class_code, "class", l_max=2500, units=None)

huble = HubbleParameter(zz, camb_code, "camb")
huble2 = HubbleParameter(zz, class_code, "class")

d_a = Angular_diameter_distance(zz, camb_code, "camb")
d_a2 = Angular_diameter_distance(zz, class_code, "class")

f1 = growth_rate(zz, camb_code, "camb")
f2 = growth_rate(zz, class_code, "class")

'''
for i in range(len(z_plot)):
    plt.loglog(p1[0], p1[2][i], label = '{}'.format(z_plot[i]))
plt.legend(title='z')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$')
plt.xlabel(r'$k\:[h/Mpc]$')
plt.show()
plt.clf()
'''
#plt.loglog(p2[0], p2[1][30])
#plt.show()

for i in range(len(z_plot)):
    plt.loglog(p2[0], p2[2][i], label = '{}'.format(z_plot[i]))

plt.legend(title='z')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$')
plt.xlabel(r'$k\:[h/Mpc]$')
plt.show()
plt.clf()

plt.plot(a1[0], a1[1])
plt.plot(a2[0], a2[1])
plt.ylabel(r'$l(l+1)\:C^{TT}_l\: / 2\pi$')
plt.xlabel('Multipole moment l')
#plt.show()
plt.clf()

plt.plot(zz,d_a)
plt.plot(zz,d_a2)
plt.ylabel(r'$d_A(z)\: [Mpc]$')
plt.xlabel('z')
#plt.show()
plt.clf()

plt.plot(zz, huble)
plt.plot(zz, huble2)
plt.ylabel(r'$H(z)\:[Km\:s⁻¹\:Mpc⁻¹]$')
plt.xlabel('z')
#plt.show()
plt.clf()

plt.plot(zz, f1)
plt.plot(zz, f2)
plt.ylabel(r'Growth rate $f_g(z)$')
plt.xlabel(r'$z$')
#plt.show()
plt.clf() 

def save_data(code, zz, hubble_array, om_b_array, om_cdm_array, APSunits=None):

    if code == "camb":
        for H in hubble_array:
            for om_b in om_b_array:
                for om_cdm in om_cdm_array:
                    camb_code = run("camb", H, om_b, om_cdm)
                    MPS = Matter_power_spectrum(zz, camb_code, "camb", l_max=2500, zplot=[0])
                    APS = Angular_power_spectrum(zz, camb_code, "camb", l_max=2500, units=None)
                    HP = HubbleParameter(zz, camb_code, "camb")
                    ADD = Angular_diameter_distance(zz, camb_code, "camb")
                    GR = growth_rate(zz, camb_code, "camb")
                    P = []
                    for i in range(len(zz)):
                        PK = MPS[1][i]
                        #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_P(k)/camb_P(k)"+ "_z="+str(np.round(zz[i],5)) +"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", PK)
                        #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_P(k)/camb_k.txt", MPS[0])
                    #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_C_l(l)/camb_Cl(l)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", APS[1])
                    #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_H(z)/camb_H(z)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", HP)
                    #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_d_A(z)/camb_d_A(z)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", ADD)
                    #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_f(z)/camb_f(z)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", GR)

    if code == "class":
            for H in hubble_array:
                for om_b in om_b_array:
                    for om_cdm in om_cdm_array:
                        class_code = run("class", H, om_b, om_cdm)
                        MPS = Matter_power_spectrum(zz, class_code, "class", l_max=2500, zplot=[0])
                        APS = Angular_power_spectrum(zz, class_code, "class", l_max=2500, units=None)
                        HP = HubbleParameter(zz, class_code, "class")
                        ADD = Angular_diameter_distance(zz, class_code, "class")
                        GR = growth_rate(zz, class_code, "class")
                        P = []
                        for i in range(len(zz)):
                            PK = MPS[1][i]
                            #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_P(k)/class_P(k)"+ "_z="+str(np.round(zz[i],5)) +"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", PK)
                            #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_P(k)/class_k.txt", MPS[0])
                        #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_C_l(l)/class_Cl(l)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", APS[1])
                        #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_H(z)/class_H(z)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", HP)
                        #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_d_A(z)/class_d_A(z)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", ADD)
                        #np.savetxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_f(z)/class_f(z)"+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", GR)
                        

'''
hubble_array = np.arange(65,75,5)
om_b_array = np.arange(0.01,0.03,0.01)
om_cdm_array = np.arange(0.1,0.3,0.1)
#save_data("camb", zz, hubble_array, om_b_array, om_cdm_array)
#save_data("class", zz, hubble_array, om_b_array, om_cdm_array)

k1 = np.loadtxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_P(k)/camb_k.txt")
k2 = np.loadtxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_P(k)/class_k.txt")

pk1 = np.loadtxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_P(k)/camb_P(k)_z=0.0_H=70_omb=0.022_omcmd=0.122.txt")
pk2 = np.loadtxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_P(k)/class_P(k)_z=0.0_H=70_omb=0.022_omcmd=0.122.txt")

h1 = np.loadtxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_H(z)/camb_H(z)_H=70_omb=0.022_omcmd=0.122.txt")
h2 = np.loadtxt("/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_H(z)/class_H(z)_H=70_omb=0.022_omcmd=0.122.txt")

plt.loglog(k1, pk1)
plt.loglog(k2, pk2)
#plt.show()

plt.plot(zz, h1)
plt.plot(zz, h2)
#plt.show()
plt.clf()
'''