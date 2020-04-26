#!/usr/bin/env python
#coding:utf8
"""
* Inputs a permittivity function of a material (by default loaded from the meep_materials.py module)
* Outputs the signal (amplitude and phase) obtained from quadratic detection of a s-SNOM microscope
(c) Filip Dominec 2013
Licensed under BSD license terms

See also:
A. Cvitkovic et al: "Analytical model for quantitative prediction of material contrasts in scattering-type 
    near-field optical microscopy", Optics Express, Vol. 15, Issue 14, pp. 8550-8565 (2007)

Huber, A.J.: "Nanoscale Surface-polariton Spectroscopy by Mid- and Far-infrared Near-field Microscopy"
    PhD thesis, 2011,  Walter Schottky Institut, Technische Universität München
    http://books.google.fr/books?id=9\_tWywAACAAJ

"""

import sys, os
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import meep_materials 

# == User settings ==
#g = .7*np.exp(.1j)
g = .95*np.exp(.1j)         ## Magic constant; value of .95*np.exp(.1j) suggested by Huber
""" observations:
        the lower amplitude of `g' shifts the resonance to lower frequencies (i.e. corresponds to additional tip-sample spacing)
        phase between -.05 and -.15 makes the curve rectangular or sharp-edged
"""
R           = 100e-9      ## tip radius
L           = 4000e-9     ## probe length (not very important for results)
dist        = 100e-9      ## the average probe distance from sample
delta_dist  = 50e-9       ## the amplitude of mechanical oscillation of the probe (delta_dist must be lower than dist)
# /User settings



def analytic_eps(mat, freq):#{{{
    """ Analytic Lorentz model (copied from meep_utils.py """
    complex_eps = mat.eps
    for polariz in mat.pol:
        complex_eps += polariz['sigma'] * polariz['omega']**2 / (polariz['omega']**2 - freq**2 - 1j*freq*polariz['gamma']) 
    return complex_eps # + sum(0)
#}}}

def finite_dipole_model(beta, g, R, H, L):#{{{
    """ 
    beta = (eps-1)/(eps+1)   .. sample "electrostatic reflection factor",
    g .. empirical coefficient,
    L .. half-length of ellipsoid,
    R .. tip radius,
    H .. tip-sample distance (along x-axis).

    """
    return beta*(g-(R+H)/L)*np.log(4*L/(4*H+3*R)) / (np.log(4*L/R) - beta*(g-(3*R+4*H)/(4*L))*np.log(2*L/(2*H+R)))
#}}}

## Dictionary of MEEP's material models
materials = {#{{{
    'TiO2':     meep_materials.material_TiO2(),
    'SiC':      meep_materials.material_SiC(),
    'SiO2':     meep_materials.material_SiO2(),
    'Au':       meep_materials.material_Au(),
    'Si':       meep_materials.material_Si_NIR(),
    'InP':      meep_materials.material_InP(),
    'GaAs':     meep_materials.material_GaAs(),
    }#}}}


if len(sys.argv)>1 and (sys.argv[1] in materials.keys()):
    material    = materials[sys.argv[1]]
    output_file     = sys.argv[1]
else:
    print "Error: the first parameter may be one of the following materials:", materials.keys()
    exit(1)

if len(sys.argv)>2 and (sys.argv[2] == "--big"):
    output_file+="_big";  plt.figure(figsize=(16,10))
else:
    plt.figure(figsize=(8,6))
freq_range=(1e12, 1e14)

freq = 10**np.arange(np.log10(freq_range[0]), np.log10(freq_range[1]), .001)
eps = analytic_eps(material, freq)
n   = (analytic_eps(material, freq)**.5).real
k   = (analytic_eps(material, freq)**.5).imag

z_raw_lo  = finite_dipole_model((eps-1)/(eps+1), g, R, 50e-9, L)
z_raw_mid  = finite_dipole_model((eps-1)/(eps+1), g, R, 100e-9, L)      ## signal at center position
z_raw_hi  = finite_dipole_model((eps-1)/(eps+1), g, R, 150e-9, L)

plt.subplot(211)
plt.title("s-SNOM signal for " + sys.argv[1])
plt.plot(freq, eps.real, color='k', label="$\\varepsilon'$")
plt.plot(freq, eps.imag, color='k', label="$\\varepsilon''$", ls='--')
plt.plot(freq, abs(z_raw_lo),     color='red',    lw=.3, label=("d-signal %d nm"  % round((dist-delta_dist)/1e-9)))
plt.plot(freq, abs(z_raw_mid),    color='green',  lw=.3, label=("d-signal %d nm"  % round((dist)/1e-9)))
plt.plot(freq, abs(z_raw_hi),    color='blue',   lw=.3, label=("d-signal %d nm"  % round((dist+delta_dist)/1e-9)))
plt.xlim(freq_range); plt.xscale('log'); plt.grid(True); plt.yscale('symlog'); plt.legend(loc=('upper left'),prop={'size':6});
plt.ylabel("Permittivity and direct signal")

plt.subplot(212)
plt.plot(freq, abs(z_raw_lo-2*z_raw_mid+z_raw_hi), color='k', marker='o', markersize=0, label="quad signal amplitude")
plt.plot(freq, np.unwrap(np.angle(z_raw_lo-2*z_raw_mid+z_raw_hi))/(2*np.pi), color='grey', marker='o', markersize=0, label="quad signal phase / $(2\\pi)$")
plt.xlim(freq_range); plt.xscale('log'); plt.grid(True); plt.legend(); plt.legend(loc=('upper left'),prop={'size':6}); 
percm = 3e10
plt.axvspan(975*percm, 1050*percm, facecolor='#FFDD00', alpha=.5, lw=0)
plt.axvspan(2.50e12, 2.54e12, facecolor='#FFDD00', alpha=1, lw=0)
plt.ylim((0,2.))
plt.ylabel("Quadratic signal")

#plt.tight_layout(True);   ## has cut the xlabel...
plt.xlabel("Frequency $f$ [Hz]\n\n");
plt.savefig("signal_"+output_file)


