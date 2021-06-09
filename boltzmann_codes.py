# -*- coding: utf-8 -*-
"""
Este script define funciones que calculan cantidades cosmologicas (matter power spectrum,
angular power spectrum, angular diameter distance, hubble parameter, growth rate) para
parametros cosmologicos a eleccion, utilizando ya sea CAMB o CLASS.
"""
import sys
reload(sys)
sys.setdefaultencoding('utf8')
import camb
from classy import Class
import numpy as np
import matplotlib.pyplot as plt

c = 3e5

def run(zz, code, Hubble0, Om_b, Om_cdm, As=2.3e-9, ns=0.96,r=0, z_pk_max=4,tau=0.09, k_max=3):
    """initialize CAMB or CLASS

    Args:
        zz (list): array of redshifts
        code (str): "camb" or "class"
        Hubble0 (float): hubble constant in [km/s/Mpc], ex: 70
        Om_b (float): Omega baryon, ex: 0.04
        Om_cdm ([type]): Omega dark matter, ex: 0.24
        As (float, optional): comoving curvature power at k=pivot_scalar. Defaults to 2.3e-9.
        ns (float, optional): scalar spectral index. Defaults to 0.96.
        r (int, optional): tensor to scalar ratio at pivot. Defaults to 0.
        z_pk_max (int, optional): maximum redshift to compute for power spectra. Defaults to 4.
        tau (float, optional): optical depth. Defaults to 0.09.
        k_max (float, optional): maximum k to calculate (not k/h). Defaults to 3.

    Returns:
        if code = camb:
            list: parameters, results
        if code = class:
            func: Class initialized with cosmology
    """
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
    """compute hubble parameter

    Args:
        z (list): array of redshifts
        run (func): initialization CAMB or CLASS
        code (str): "camb" or "class"

    Returns:
        list: hubble parameter
    """
    if code == "camb":
        pars, results = run[0], run[1]
        A = results.get_BAO(z, pars)
        H = A[:,1] 
        return H

    if code == "class":
        hubble_vect = np.vectorize(run.Hubble)
        H = hubble_vect(z)*c
        return H

def Angular_diameter_distance(z, run, code):
    """compute angular diameter distance

    Args:
        z (list): array of redshifts
        run (func): initialization CAMB or CLASS
        code (str): "camb" or "class"

    Returns:
        list: angular diameter distance
    """
    if code == "camb":
        pars, results = run[0], run[1]
        A = results.get_BAO(z, pars)
        D_A = A[:,2]
        return D_A

    if code == "class":
        DA_vect = np.vectorize(run.angular_distance)
        D_A = DA_vect(z)
        return D_A

def growth_rate(z, run, code, min_k=1e-4, max_k=3, n_points = 2000, cmb_units='muK'):
    """compute growth rate

    Args:
        z (list): array of redshifts
        run (func): initialization CAMB or CLASS
        code (str): "camb" or "class"
        min_k (float, optional): minimum k to calculate. Defaults to 1e-4.
        max_k (float, optional): maximum k to calculate. Defaults to 3.
        n_points (int, optional): lenght of k array. Defaults to 2000.

    Returns:
        list: growth rate
    """
    if code == "camb":
        pars, results = run[0], run[1]
        results.get_cmb_power_spectra(pars, CMB_unit=cmb_units)
        results.get_matter_power_spectrum(minkh=min_k, maxkh=max_k, npoints = n_points)
        fsigma8, sigma8 = results.get_fsigma8(), results.get_sigma8()
        f_growth_rate = np.flip(fsigma8/sigma8)
        return f_growth_rate

    if code == "class":
        f_vect = np.vectorize(run.scale_independent_growth_factor_f)
        growth_rate_f = f_vect(z)
        return growth_rate_f

def Angular_power_spectrum(z, run, code, l_max, units, lens_accuracy=0):
    """compute angular power spectrum

    Args:
        z (list): array of redshifts
        run (func): initialization CAMB or CLASS
        code (str): "camb" or "class"
        l_max (float): maximum multipolar moment
        units (str): None or 'muK'
        lens_accuracy (float, optional): Set to 1 or higher if you want to get the lensing potential accurate. Defaults to 0.

    Returns:
        list: angular power spectrum
    """
    if code == "camb":
        pars, results = run[0], run[1]
        pars.set_for_lmax(l_max, lens_potential_accuracy=lens_accuracy)
        powers =results.get_cmb_power_spectra(pars, CMB_unit= units)
        totCl=powers['total']
        ls = np.arange(totCl.shape[0])
        Cl_TT = totCl[:,0]
        #Cl_EE = totCl[:,1]
        #Cl_BB = totCl[:,2]
        #Cl_ET = totCl[:,3]
        return ls, Cl_TT

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
        return ls, Cl_TT


def P_k_camb(z, pk, zplot, P_z, units,h):
    """compute matter power spectrum for zplot redshifts
       function used in func Matter_power_spectrum

    Args:
        z (list): array of redshifts
        pk (list): matter power spectrum array computed by camb o class
        zplot (list): redshifts to return
        P_z (list): Pk(zplot)
        units (str): 'Mpc/h' or 'Mpc'
        h (float): hubble constant
    """
    if units == '(Mpc/h)^3':
        for i in range(0, len(z)):
            for z_j in zplot:
                if z[i] == z_j:
                    P_z.append(pk[i,:])
    if units == 'Mpc^3':
        for i in range(0, len(z)):
            for z_j in zplot:
                if z[i] == z_j:
                    P_z.append(pk[i,:])

def Matter_power_spectrum(z, run, code, l_max, zplot, kmin=1e-4, kmax=2, number_points=2000, plot="yes", units='(Mpc/h)^3'):
    """compute matter power spectrum

    Args:
        z (list): array of redshifts
        run (func): initialization CAMB or CLASS
        code (str): "camb" or "class"
        l_max (float): maximum multipolar moment
        zplot (list): redshifts to return
        kmin (float, optional): minimum k to calculate. Defaults to 1e-4.
        kmax (float, optional): maximum k to calculate. Defaults to 2.
        number_points (int, optional): lenght of k array. Defaults to 2000.
        plot (str, optional): if 'yes' compute P_k for zplot redshifts
                              if 'no' compute P_k for z redshifts. Defaults to "yes".
        units (str, optional): '(Mpc/h)^3' or 'Mpc^3'. Defaults to '(Mpc/h)^3'.

    Returns:
        list: k, P_k
    """
    if code == "camb":
        if units == '(Mpc/h)^3':
            pars = run[0]
            results = camb.get_results(pars)
            h = results.h_of_z(0)
            kh,z, pk = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = number_points)
            P_z = []
            if plot == "yes":
                P_k_camb(z, pk, zplot, P_z, units,h)
                return kh, pk, P_z
            if plot == "no":
                return kh, pk
    
        if units == 'Mpc^3':
            pars = run[0]
            results = camb.get_results(pars)
            h = results.h_of_z(0)
            #kh,z, pk = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints = number_points)
            kh,z,pk = results.get_linear_matter_power_spectrum(hubble_units=False, k_hunit=False) 
            P_z = []
            if plot == "yes":
                P_k_camb(z, pk, zplot, P_z, units,h)
                return kh, pk, P_z
            if plot == "no":
                return kh, pk


    if code == "class":
        if units == '(Mpc/h)^3':
            kk = np.linspace(kmin, kmax, number_points)
            Pk = []
            Pk_z = []
            h = run.h()
            pk_vect = np.vectorize(run.pk)
            for zz in z:
                P = pk_vect(kk*h,zz) * h**3
                Pk.append(P)
            return kk, Pk

        if units == 'Mpc^3':
            h = run.h()
            kmin = kmin/h
            kmax = kmax/h
            kk = np.linspace(kmin, kmax, number_points)
            Pk = []
            Pk_z = []
            k = kk
            pk_vect = np.vectorize(run.pk)
            for zz in z:
                P = pk_vect(kk,zz)
                Pk.append(P)
            return k, Pk

def params_eps(save_data, params, epsilon, code, zz, lmax, z_pk=[0], perturbations=True, APSunits=None, MPSunits='(Mpc/h)^3',
               z_pk_max=4, k_max=3, lens_accuracy=0):
    for i in range(len(params)):
        if params[i] == 0:
            params[i] = epsilon
        params[i] = params[i] * (1+epsilon)
        if i>0:
            params[i-1] = params[i-1] / (1+epsilon)
        data = save_data(code, zz, params[0], params[1], params[2], lmax, z_pk, perturbations,
                  APSunits, MPSunits, params[3], params[4], params[5], z_pk_max, 
                  params[6], k_max, lens_accuracy)
        print(params)
    params[-1] = params[-1] / (1+epsilon)
    params_plus_eps = np.array(params) * (1+epsilon)
        #print(data)
    #print(params)
    #print(params_plus_eps)
    return params_plus_eps, data

def save_data(code, zz, Hubble0, om_b, om_cdm, lmax, z_pk=[0], perturbations=True, APSunits=None, MPSunits='(Mpc/h)^3',
              As=2.3e-9, ns=0.96,r=0, z_pk_max=4,tau=0.09, k_max=3, lens_accuracy=0):

        if code == "camb":
            camb_code = run(zz, "camb", Hubble0, om_b, om_cdm, As=As, ns=ns, z_pk_max=z_pk_max, tau=tau, k_max=k_max)
            h = HubbleParameter(zz, camb_code, "camb")
            da = Angular_diameter_distance(zz, camb_code, "camb")
            if perturbations == True:
                f = growth_rate(zz, camb_code, "camb", cmb_units=APSunits)
                pk = Matter_power_spectrum(zz, camb_code, "camb", l_max=lmax, zplot=z_pk, units= MPSunits, plot='yes')
                cl = Angular_power_spectrum(zz, camb_code, "camb", l_max=lmax, units=APSunits, lens_accuracy=lens_accuracy)
                k = pk[0]
                ls = cl[0]
                return h, da, f, [k, pk[2]], [ls, cl]
            if perturbations != True:
                return h, da

        if code == "class":
            class_code = run(zz, "class", Hubble0, om_b, om_cdm)
            h = HubbleParameter(zz, class_code, "class")
            da = Angular_diameter_distance(zz, class_code, "class")
            if perturbations == True:
                f = growth_rate(zz, class_code, "class") 
                pk = Matter_power_spectrum(zz, class_code, "class", l_max=lmax, zplot=z_pk, units=MPSunits, plot='yes')
                cl = Angular_power_spectrum(zz, class_code, "class", l_max=lmax, units=APSunits) 
                k = pk[0]
                ls = cl[0]
                return  h, da, f, [k, pk[2]], [ls, cl]
            if perturbations != True:
                return h, da
                

def save_data_array(code, zz, hubble_array, om_b_array, om_cdm_array, lmax, perturbations='True', APSunits=None, MPSunits='(Mpc/h)^3',
              As=2.3e-9, ns=0.96,r=0, z_pk_max=4,tau=0.09, k_max=3, lens_accuracy=0):

    """Calculate cosmological quantities for all combinations of hubble constant,
       omega baryion and omega dark matter.

    Args:
        code (str): "camb" or "class"
        zz (list): array of redshifts
        hubble_array (list): array of hubble constants
        om_b_array (list): array of omega baryion
        om_cdm_array (list): array of omega dark matter to calculate
        l_max (float): maximum multipolar moment
        perturbations (str): 'True' to compute Pk and Cl, 'False' only compute H, d_A, f
        APSunits (str, optional): angular power spectrum units. None or "muK". Defaults to None.
        MPSunits (str, optional): matter power spectrum units. '(Mpc/h)^3' o 'Mpc^3'. Defaults to '(Mpc/h)^3'.
        As (float, optional): comoving curvature power at k=pivot_scalar. Defaults to 2.3e-9.
        ns (float, optional): scalar spectral index. Defaults to 0.96.
        r (int, optional): tensor to scalar ratio at pivot. Defaults to 0.
        z_pk_max (int, optional): maximum redshift to compute for power spectra. Defaults to 4.
        tau (float, optional): optical depth. Defaults to 0.09.
        k_max (float, optional): maximum k to calculate (not k/h). Defaults to 3.
        lens_accuracy (float, optional): Set to 1 or higher if you want to get the lensing potential accurate. Defaults to 0.

    Returns:
        list: array of [k,Pk], [l,Cl], H, d_A, f
    """
    MPS = []
    APS = []
    HP = []
    ADD = []
    GR = []
    k = []

    if code == "camb":
        for H in hubble_array:
            for om_b in om_b_array:
                for om_cdm in om_cdm_array:
                    camb_code = run(zz, "camb", H, om_b, om_cdm, As=As, ns=ns, z_pk_max=z_pk_max, tau=tau, k_max=k_max)
                    pk = Matter_power_spectrum(zz, camb_code, "camb", l_max=lmax, zplot=[0], units= MPSunits)
                    cl = Angular_power_spectrum(zz, camb_code, "camb", l_max=lmax, units=APSunits, lens_accuracy=lens_accuracy)
                    h = HubbleParameter(zz, camb_code, "camb")
                    da = Angular_diameter_distance(zz, camb_code, "camb")
                    f = growth_rate(zz, camb_code, "camb")
                    P_z = []
                    for i in range(len(zz)):
                        P = pk[1][i]
                        P_z.append(P)
                    MPS.append(P_z)
                    APS.append(cl[1])
                    HP.append(h)
                    ADD.append(da)
                    GR.append(f)
                    k = pk[0]
                    ls = cl[0]
        MPS = [k, MPS]
        APS = [ls, APS]

        return MPS, APS, HP, ADD, GR

    if code == "class":
        for H in hubble_array:
            for om_b in om_b_array:
                for om_cdm in om_cdm_array:
                    class_code = run(zz, "class", H, om_b, om_cdm)
                    pk = Matter_power_spectrum(zz, class_code, "class", l_max=lmax, zplot=[0], units=MPSunits)
                    cl = Angular_power_spectrum(zz, class_code, "class", l_max=lmax, units=APSunits)
                    h = HubbleParameter(zz, class_code, "class")
                    da = Angular_diameter_distance(zz, class_code, "class")
                    f = growth_rate(zz, class_code, "class")
                    P_z = []
                    for i in range(len(zz)):
                        P = pk[1][i]
                        P_z.append(P)
                    MPS.append(P_z)
                    APS.append(cl[1])
                    HP.append(h)
                    ADD.append(da)
                    GR.append(f)
                    k = pk[0]
                    ls = cl[0]
        MPS = [k, MPS]
        APS = [ls, APS]

        return MPS, APS, HP, ADD, GR

def save_data_txt(code, zz, hubble_array, om_b_array, om_cdm_array, paths, names,
                  lmax,APSunits=None, MPSunits='(Mpc/h)^3',
                  As=2.3e-9, ns=0.96,r=0, z_pk_max=4,tau=0.09, k_max=3, lens_accuracy=0):
    """Save cosmological quantities for all combinations of hubble constant,
       omega baryion and omega dark matter in .txt.

    Args:
        code (str): "camb" or "class"
        zz (list): array of redshifts
        hubble_array (list): array of hubble constants
        om_b_array (list): array of omega baryion
        om_cdm_array (list): array of omega dark matter to calculate
        paths (list[str]): array of paths to save .txt. order: PK,k,Cl,H,d_A,f,l
        names (lisrt[str]): array of names of .txt. order: PK,k,Cl,H,d_A,f,l
        l_max (float): maximum multipolar moment
        APSunits (str, optional): angular power spectrum units. None or "muK". Defaults to None.
        MPSunits (str, optional): matter power spectrum units. '(Mpc/h)^3' o 'Mpc^3'. Defaults to '(Mpc/h)^3'.
        As (float, optional): comoving curvature power at k=pivot_scalar. Defaults to 2.3e-9.
        ns (float, optional): scalar spectral index. Defaults to 0.96.
        r (int, optional): tensor to scalar ratio at pivot. Defaults to 0.
        z_pk_max (int, optional): maximum redshift to compute for power spectra. Defaults to 4.
        tau (float, optional): optical depth. Defaults to 0.09.
        k_max (float, optional): maximum k to calculate (not k/h). Defaults to 3.
        lens_accuracy (float, optional): Set to 1 or higher if you want to get the lensing potential accurate. Defaults to 0.
    """

    if code == "camb":
        for H in hubble_array:
            for om_b in om_b_array:
                for om_cdm in om_cdm_array:
                    camb_code = run(zz, "camb", H, om_b, om_cdm)
                    MPS = Matter_power_spectrum(zz, camb_code, "camb", l_max=2500, zplot=[0])
                    APS = Angular_power_spectrum(zz, camb_code, "camb", l_max=2500, units=None)
                    HP = HubbleParameter(zz, camb_code, "camb")
                    ADD = Angular_diameter_distance(zz, camb_code, "camb")
                    GR = growth_rate(zz, camb_code, "camb")
                    P = []
                    for i in range(len(zz)):
                        PK = MPS[1][i]
                        np.savetxt(paths[0][0] + names[0][0] + "_z="+str(np.round(zz[i],5)) +"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", PK)
                        np.savetxt(paths[0][1]+ names[0][1]+ ".txt", MPS[0])                        
                    np.savetxt(paths[0][2] + names[0][2]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", APS[1])
                    np.savetxt(paths[0][3] + names[0][3]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", HP)
                    np.savetxt(paths[0][4] + names[0][4]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", ADD)
                    np.savetxt(paths[0][5] + names[0][5]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", GR)
                    np.savetxt(paths[0][6]+ names[0][6]+ ".txt", APS[0])

    if code == "class":
            for H in hubble_array:
                for om_b in om_b_array:
                    for om_cdm in om_cdm_array:
                        class_code = run(zz, "class", H, om_b, om_cdm)
                        MPS = Matter_power_spectrum(zz, class_code, "class", l_max=2500, zplot=[0])
                        APS = Angular_power_spectrum(zz, class_code, "class", l_max=2500, units=None)
                        HP = HubbleParameter(zz, class_code, "class")
                        ADD = Angular_diameter_distance(zz, class_code, "class")
                        GR = growth_rate(zz, class_code, "class")
                        P = []
                        for i in range(len(zz)):
                            PK = MPS[1][i]
                            np.savetxt(paths[1][0] + names[1][0]+ "_z="+str(np.round(zz[i],5)) +"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", PK)
                            np.savetxt(paths[1][1] + names[1][1]+".txt", MPS[0])
                        np.savetxt(paths[1][2] + names[1][2]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", APS[1])
                        np.savetxt(paths[1][3]+ names[1][3]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", HP)
                        np.savetxt(paths[1][4]+ names[1][4]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", ADD)
                        np.savetxt(paths[1][5]+ names[1][5]+"_H="+str(H)+"_omb="+str(om_b)+"_omcmd="+str(om_cdm)+ ".txt", GR)
                        np.savetxt(paths[1][6]+ names[1][6]+ ".txt", APS[0])
'''
z_plot = [0,1,2]
zmin = 0
zmax = 3
step = 0.05
zz = np.arange(zmin,zmax,step)
#camb_code = run("camb", 70, 0.022, 0.122)
#class_code = run("class", 70, 0.022, 0.122)

camb_code = run(zz, "camb", 70, 0.04, 0.24)
class_code = run(zz, "class", 70, 0.04, 0.24)

p1 = Matter_power_spectrum(z_plot, camb_code, "camb", l_max=2500, zplot=z_plot, plot='yes', units='Mpc^3')
p2 = Matter_power_spectrum(z_plot, class_code, "class", l_max=2500, zplot=z_plot, plot='yes', units='Mpc^3')

a1 = Angular_power_spectrum(zz, camb_code, "camb", l_max=2500, units=None)
a2 = Angular_power_spectrum(zz, class_code, "class", l_max=2500, units=None)

huble = HubbleParameter(zz, camb_code, "camb")
huble2 = HubbleParameter(zz, class_code, "class")

d_a = Angular_diameter_distance(zz, camb_code, "camb")
d_a2 = Angular_diameter_distance(zz, class_code, "class")

f1 = growth_rate(zz, camb_code, "camb")
f2 = growth_rate(zz, class_code, "class")


for i in range(len(z_plot)):
    plt.loglog(p1[0], p1[2][i], label = '{}'.format(z_plot[i]))
plt.legend(title='z')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$')
plt.xlabel(r'$k\:[h/Mpc]$')
plt.show()
plt.clf()

#plt.loglog(p2[0], p2[1][30])
#plt.show()

for i in range(len(z_plot)):
    plt.loglog(p2[0], p2[1][i], label = '{}'.format(z_plot[i]))

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