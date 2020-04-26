#!/usr/bin/env python
#coding:utf8
import numpy
def finite_dipole_model(beta, g, R, H, L):#{{{
    """ 
    Calculates the scattering efficiency of a s-SNOM microscope, see also:
    A. Cvitkovic et al: "Analytical model for quantitative prediction of material contrasts in scattering-type 
        near-field optical microscopy",
    Optics Express, Vol. 15, Issue 14, pp. 8550-8565 (2007)

    beta = (eps-1)/(eps+1)   .. sample "electrostatic reflection factor",
    g .. empirical coefficient,
    L .. half-length of ellipsoid,
    R .. tip radius,
    H .. tip-sample distance (along x-axis).

    """
    return beta*(g-(R+H)/L)*numpy.log(4*L/(4*H+3*R)) / (numpy.log(4*L/R) - beta*(g-(3*R+4*H)/(4*L))*numpy.log(2*L/(2*H+R)))
#}}}
