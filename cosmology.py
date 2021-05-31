# -*- coding: utf-8 -*-
"""
Este script utiliza funciones, importadas de boltzmann_codes.py, que calculan cantidades 
cosmologicas (matter power spectrum, angular power spectrum, angular diameter distance, 
hubble parameter, growth rate) para parametros cosmologicos a eleccion, utilizando ya sea 
CAMB o CLASS
"""
import sys
reload(sys)
sys.setdefaultencoding('utf8')
import numpy as np
import matplotlib.pyplot as plt
path = '/home/gerald/Documentos/proyecto_cmb/cmb-1/boltzmann_codes.py' #boltzmann_codes.py path
sys.path.insert(0, path)
import boltzmann_codes

z_plot = [0,1,2]
zmin = 0
zmax = 3
step = 0.05
zz = np.arange(zmin,zmax,step)
camb_code = boltzmann_codes.run(zz, "camb", 70, 0.04, 0.24)
class_code = boltzmann_codes.run(zz, "class", 70, 0.04, 0.24)

p1 = boltzmann_codes.Matter_power_spectrum(z_plot, camb_code, "camb", l_max=2500, zplot=z_plot, plot='yes', units='(Mpc/h)^3')
p2 = boltzmann_codes.Matter_power_spectrum(z_plot, class_code, "class", l_max=2500, zplot=z_plot, plot='yes', units='(Mpc/h)^3')

a1 = boltzmann_codes.Angular_power_spectrum(zz, camb_code, "camb", l_max=2500, units=None)
a2 = boltzmann_codes.Angular_power_spectrum(zz, class_code, "class", l_max=2500, units=None)

huble = boltzmann_codes.HubbleParameter(zz, camb_code, "camb")
huble2 = boltzmann_codes.HubbleParameter(zz, class_code, "class")

d_a = boltzmann_codes.Angular_diameter_distance(zz, camb_code, "camb")
d_a2 = boltzmann_codes.Angular_diameter_distance(zz, class_code, "class")

f1 = boltzmann_codes.growth_rate(zz, camb_code, "camb")
f2 = boltzmann_codes.growth_rate(zz, class_code, "class")

hubble_array = np.arange(65,75,5)
om_b_array = np.arange(0.03,0.05,0.01)
om_cdm_array = np.arange(0.23,0.25,0.01)

'''
paths = [["/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_P(k)/","/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_P(k)/", 
"/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_C_l(l)/","/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_H(z)/",
"/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_d_A(z)/", "/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_f(z)/",
"/home/gerald/Documentos/proyecto_cmb/cmb-1/data/camb/camb_C_l(l)/"],
["/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_P(k)/", "/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_P(k)/",
"/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_C_l(l)/", "/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_H(z)/",
"/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_d_A(z)/", "/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_f(z)/",
"/home/gerald/Documentos/proyecto_cmb/cmb-1/data/class/class_C_l(l)/"]]

names = [["camb_P(k)", "camb_k", "camb_Cl(l)", "camb_H(z)", "camb_d_A(z)", "camb_f(z)", "camb_ls"],
["class_P(k)", "class_k", "class_Cl(l)", "class_H(z)", "class_d_A(z)", "class_f(z)", "class_ls"]]

boltzmann_codes.save_data_txt("camb", zz, hubble_array, om_b_array, om_cdm_array, paths, names, lmax=2500)
boltzmann_codes.save_data_txt("class", zz, hubble_array, om_b_array, om_cdm_array, paths, names,lmax=2500)
'''

#save_camb = boltzmann_codes.save_data("camb", zz, hubble_array, om_b_array, om_cdm_array,lmax=2500)
#save_class = boltzmann_codes.save_data("class", zz, hubble_array, om_b_array, om_cdm_array,lmax=2500)

for i in range(len(z_plot)):
    plt.loglog(p1[0], p1[2][i], label = '{}'.format(z_plot[i]))
plt.legend(title='z')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$')
plt.xlabel(r'$k\:[h/Mpc]$')
plt.show()
plt.clf()

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
plt.show()
plt.clf()

plt.plot(zz,d_a)
plt.plot(zz,d_a2)
plt.ylabel(r'$d_A(z)\: [Mpc]$')
plt.xlabel('z')
plt.show()
plt.clf()

plt.plot(zz, huble)
plt.plot(zz, huble2)
plt.ylabel(r'$H(z)\:[Km\:s⁻¹\:Mpc⁻¹]$')
plt.xlabel('z')
plt.show()
plt.clf()

plt.plot(zz, f1)
plt.plot(zz, f2)
plt.ylabel(r'Growth rate $f_g(z)$')
plt.xlabel(r'$z$')
plt.show()
plt.clf() 