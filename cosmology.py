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
import time

t0= time.clock()
#print(t0)

z_plot = [0,1]
zmin = 0
zmax = 3
step = 0.05
zz = np.arange(zmin,zmax,step)

camb_code = boltzmann_codes.run(zz, "camb", 70, 0.04, 0.24)
class_code = boltzmann_codes.run(zz, "class", 70, 0.04, 0.24)

#boltzmann_codes.params_eps(boltzmann_codes.save_data, [[70,False], [0.04,False], [0.24,False], [2.3e-9, False],
#                           [0.96,True], [0,False], [0.09,False]], 0.01, 'camb', zz, lmax=2500, perturbations=True, z_pk=zz)


#hubble_array = np.arange(65,75,7.5)
#om_b_array = np.arange(0.03,0.05,0.015)
#om_cdm_array = np.arange(0.23,0.25,0.015)
#save_camb_array = boltzmann_codes.save_data("camb", zz, 70, 0.04, 0.24, lmax=2500, perturbations=False)
#save_class_array = boltzmann_codes.save_data("class", zz, 70, 0.04, 0.24, lmax=2500, perturbations=False)
#save_class = boltzmann_codes.save_data_array("class", zz, hubble_array, om_b_array, om_cdm_array,lmax=2500)

t = time.clock() - t0
#print(t)


p1 = boltzmann_codes.Matter_power_spectrum(zz, camb_code, "camb", l_max=2500, zplot=z_plot, plot='yes', units='(Mpc/h)^3')
p2 = boltzmann_codes.Matter_power_spectrum(zz, class_code, "class", l_max=2500, zplot=z_plot, plot='yes', units='(Mpc/h)^3')


a1 = boltzmann_codes.Angular_power_spectrum(zz, camb_code, "camb", l_max=2500, units=None)
a2 = boltzmann_codes.Angular_power_spectrum(zz, class_code, "class", l_max=2500, units=None)


huble = boltzmann_codes.HubbleParameter(zz, camb_code, "camb")
huble2 = boltzmann_codes.HubbleParameter(zz, class_code, "class")


d_a = boltzmann_codes.Angular_diameter_distance(zz, camb_code, "camb")
d_a2 = boltzmann_codes.Angular_diameter_distance(zz, class_code, "class")


f1 = boltzmann_codes.growth_rate(zz, camb_code, "camb")
f2 = boltzmann_codes.growth_rate(zz, class_code, "class")

ls = iter(['solid', 'dashed', 'dashdot'])
ls2 = iter(['solid', 'dashed', 'dashdot', 'dashdot'])
error = np.ones(len(zz))*0.005
error2 = np.ones(len(p1[0]))*0.05

for i in range(len(z_plot)-1):
    plt.plot(p1[0], (p1[2][i] - p2[1][i])/p1[2][i] , label = 'z = {}'.format(z_plot[i]), color='darkviolet')
    print((p1[2][i] - p2[1][i])/p1[2][i])
plt.plot(p1[0], error2, color='k', ls='dashed')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$', fontsize=13)
plt.xlabel(r'$k\:[h/Mpc]$', fontsize=13)
plt.xlim((0.001,2))
plt.ylim((0,1.1))
#plt.show()
plt.clf()

'''plt.plot(a1[0][3:], (a1[1][3:] - a2[1][3:])/a1[1][3:])
plt.ylabel(r'$l(l+1)\:C^{TT}_l\: / 2\pi$', fontsize=13)
plt.xlabel('Multipole moment l', fontsize=13)
plt.legend()
plt.show()
plt.clf()'''

plt.plot(zz,np.abs((d_a-d_a2)/d_a), color='darkviolet')
#plt.plot(zz, error, color='k', ls='dashed')
plt.ylabel(r'$d_A(z)\: [Mpc]$', fontsize=13)
plt.xlabel(r'$z$', fontsize=13)
plt.legend()
#plt.ylim((-0.01,0.01))
#plt.show()
plt.clf()

plt.plot(zz, np.abs((huble-huble2)/huble), color='darkviolet')
#plt.plot(zz, error, color='k', ls='dashed')
plt.ylabel(r'$H(z)\:[Km\:s⁻¹\:Mpc⁻¹]$', fontsize=13)
plt.xlabel(r'$z$',fontsize=13)
plt.legend()
#plt.ylim((-0.01,0.01))
#plt.show()
plt.clf()

plt.plot(zz, np.abs((f1-f2)/f1), color='darkviolet')
#plt.plot(zz, error, color='k', ls='dashed')
plt.ylabel(r'Growth rate $f_g(z)$',fontsize=13)
plt.xlabel(r'$z$', fontsize=13)
plt.legend()
#plt.ylim((-0.01,0.01))
#plt.show()
plt.clf()

fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(zz, np.abs((huble-huble2)/huble), color='darkviolet')
axs[0, 0].set_xlabel(r'$z$', fontsize=11)
#axs[0, 0].set_ylabel(r'$H(z)$', fontsize=11)
axs[0, 0].ticklabel_format(axis='both', style='sci', scilimits=(0,0))
axs[0, 0].set_title(r'$H(z)$')
axs[1, 0].plot(zz, np.abs((d_a-d_a2)/d_a), color='darkviolet')
axs[1, 0].set_xlabel(r'$z$', fontsize=11)
#axs[1, 0].set_ylabel(r'$d_A(z)$', fontsize=11)
axs[1, 0].ticklabel_format(axis='both', style='sci', scilimits=(0,0))
axs[1, 0].set_title(r'$d_A(z)$')
axs[0, 1].plot(zz, np.abs((f1-f2)/f1), color='darkviolet')
axs[0, 1].set_xlabel(r'$z$', fontsize=11)
#axs[0, 1].set_ylabel(r'$f_g(z)$',fontsize=11)
axs[0, 1].ticklabel_format(axis='both', style='sci', scilimits=(0,0))
axs[0, 1].set_title(r'$f_g(z)$')
for i in range(len(z_plot)-1):
    axs[1, 1].plot(p1[0], (p1[2][i] - p2[1][i])/p1[2][i] , label = 'z = {}'.format(z_plot[i]), color='darkviolet')
axs[1, 1].plot(p1[0], error2, color='k', ls='dashed')
axs[1, 1].set_xlabel(r'$k\:[h/Mpc]$', fontsize=11)
#axs[1, 1].set_ylabel(r'$P(k)$', fontsize=11)
axs[1, 1].set_xlim((0.01,2))
axs[1, 1].set_ylim((0,1.1))
axs[1, 1].ticklabel_format(axis='both', style='sci', scilimits=(0,0))
axs[1, 1].set_title(r'$P(k)$')
fig.tight_layout(w_pad=4)
fig.suptitle('Relative Errors', y=0.99, fontsize=13)
plt.subplots_adjust(left=0.1, right=0.95, bottom=0.1, top=0.9)
plt.savefig('new_plots/relative_errors')
plt.show()


'''
ls = iter(['solid', 'dashed', 'dashdot'])
ls2 = iter(['solid', 'dashed', 'dashdot', 'dashdot'])

for i in range(len(z_plot)-1):
    plt.loglog(p1[0], p1[2][i], label = 'z = {}'.format(z_plot[i]), color='darkblue', ls=next(ls))
#plt.legend(title='z')
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$', fontsize=13)
plt.xlabel(r'$k\:[h/Mpc]$', fontsize=13)
#plt.show()
#plt.clf()


for i in range(len(z_plot)-1):
    plt.loglog(p2[0], p2[1][i], label = 'z = {}'.format(z_plot[i]), color='red', ls=next(ls2))

plt.legend(title='CAMB             CLASS', ncol=2)
plt.ylabel(r'$P(k)\:[(Mpc/h)³]$', fontsize=13)
plt.xlabel(r'$k\:[h/Mpc]$', fontsize=13)
plt.savefig('new_plots/P(k)')
plt.show()
plt.clf()

plt.plot(a1[0][3:], a1[1][3:], color='darkblue', label='CAMB')
plt.plot(a2[0][3:], a2[1][3:], color='red', label='CLASS')
plt.ylabel(r'$l(l+1)\:C^{TT}_l\: / 2\pi$', fontsize=13)
plt.xlabel('Multipole moment l', fontsize=13)
plt.legend()
plt.savefig('new_plots/C(l)')
plt.show()
plt.clf()

plt.plot(zz,d_a, color='darkblue', label='CAMB')
plt.plot(zz,d_a2, color='red', label='CLASS')
plt.ylabel(r'$d_A(z)\: [Mpc]$', fontsize=13)
plt.xlabel(r'$z$', fontsize=13)
plt.legend()
plt.savefig('new_plots/d_A')
plt.show()
plt.clf()

plt.plot(zz, huble, color='darkblue', label='CAMB')
plt.plot(zz, huble2, color='red', label='CLASS')
plt.ylabel(r'$H(z)\:[Km\:s⁻¹\:Mpc⁻¹]$', fontsize=13)
plt.xlabel(r'$z$',fontsize=13)
plt.legend()
plt.savefig('new_plots/H(z)')
plt.show()
plt.clf()

plt.plot(zz, f1, color='darkblue', label='CAMB')
plt.plot(zz, f2, color='red', label='CLASS')
plt.ylabel(r'Growth rate $f_g(z)$',fontsize=13)
plt.xlabel(r'$z$', fontsize=13)
plt.legend()
plt.savefig('new_plots/f(z)')
plt.show()
plt.clf() 
'''

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
'''