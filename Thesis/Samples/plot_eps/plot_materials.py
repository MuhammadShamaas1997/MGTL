#!/usr/bin/env python
#coding:utf8

import os
import sys
import traceback
import numpy as np
import matplotlib
import matplotlib.cm
import matplotlib.pyplot as plt
import meep_materials 


## Interesting materials: 
##                  was_plotted         has_model           See_also
##      Si          OK                  (poor)
##      SiC		    OK                  OK
##      SiO2		narrow              OK                  Huber?
##      TiO2		narrow              OK                  Baumard77, ..
##      InP         OK                  (questionable)
##      GaAs        OK-lim
##      InSb        OK-lim
##      Au		    OK                  OK

##      FeS2        N/A
##      Si-doped                                            Si plasmons: http://proj.ncku.edu.tw/research/articles/e/20090828/2.html

##      SrTiO3      TODO


## Analytic Lorentz model (copied from meep_utils.py)
def analytic_eps(mat, freq):#{{{
    complex_eps = mat.eps
    for polariz in mat.pol:
        complex_eps += polariz['sigma'] * polariz['omega']**2 / (polariz['omega']**2 - freq**2 - 1j*freq*polariz['gamma']) 
    return complex_eps # + sum(0)
#}}}

## Functions loading data from different format files 
## Note: you must obtain the SOPRA and SCOUT databases from the web, if you want to use them
def load_SCOUT_permittivity(filename):#{{{
    """ Reads the permittivity function from a given file with SCOUT binary format

    The SCOUT database of materials is supplied with the SCOUT program and may be freely downloaded
    from http://www.mtheiss.com/download/scout.zip
    
    Different files use different units for the x-axis (179 files with "eV", 118 "micron", 49 "1/cm"), so
    we have to build the frequency axis accordingly. The data were randomly verified against the SOPRA and luxpop data.
    """


    ## Open the file for binary access
    f = open(filename, "rb")

    ## Load the number of data points, type of x axis and its boundaries
    f.seek(151);    datalength = np.fromfile(f, dtype=np.uint16, count=1)[0]
    print(datalength)
    f.seek(160);    x_axis_type =  np.fromfile(f, dtype=np.uint8, count=1)[0]
    print(x_axis_type)
    f.seek(166);    x_start, x_end = np.fromfile(f, dtype=np.float32, count=2)
    print(x_start)
    print(x_end)
    ## Load the n, k data
    f.seek(174);    raw_eps = np.fromfile(f, dtype=np.float32, count=datalength*2)
    f.close

    eps = raw_eps[::2] + 1j*raw_eps[1::2]
 
    from scipy.constants import h, c, eV
    if x_axis_type == 2:            # 'eV'            
        freq = np.linspace(x_start*eV/h, x_end*eV/h, datalength)   
    elif x_axis_type == 3:          # 'um'           
        wavelength = np.linspace(x_start*1e-6, x_end*1e-6, datalength)
        freq = c/wavelength
    elif x_axis_type == 0:          # 'cm-1'      
        freq = np.linspace(x_start*100*c, x_end*100*c, datalength)   

    return freq, eps
#}}}
def load_SOPRA_permittivity(filename):#{{{
    data = []
    with open(filename) as f:
        for line in f.readlines():
            if line[:5] == 'DATA1': data.append(map(lambda x: float(x), line.split('*')[2:5]))
    wl, n, k = np.array(data).T
    eps = ((n+1j*k)**2)[::-1]
    freq = (2.998e8 / (wl*1e-9))[::-1]
    return freq, eps
#}}}
def load_n_k_eps(filename):#{{{
    lambda_angstr, n, k = np.loadtxt(data_source, usecols=[0,1,2], unpack=True, comments=';')
    eps = (n+1j*k)**2
    freq = 2.997e8 / (lambda_angstr*1e-10)
    return freq, n, k, eps
#}}}
def load_eps(filename):#{{{
    return freq, eps
#}}}

## List of sources for given materials
## == TiO2  ==#{{{
TiO2_files = [
        meep_materials.material_TiO2(where = None),
        meep_materials.material_TiO2(where = None, extraordinary=1.),      ## (rutile is anisotropic)
        meep_materials.material_TiO2(where = None, extraordinary=0.),
        "sopra/TIO2.MAT",
        "sopra/TIO2B.MAT",
        "luxpop/TiO2.nk",
        "luxpop/TiO2-e.nk",
        "luxpop/TiO2_llnl_cxro.nk",
        "luxpop/TiO2-e_palik.nk",
        "scout/TiO2 (amorph).b",
        "scout/TiO2 I (Jellison).b",
        "scout/TiO2 II (Jellison).b",
        "scout/TiO2 II (Jellison).b",
        "other/TiO2-Mounaix_polycryst.nk",
        ]
#}}}
## == SiC == #{{{
SiC_files = [
        meep_materials.material_SiC(where = None),
        "sopra/SIC.MAT", 
        "luxpop/SiC_llnl_cxro.nk",
        "luxpop/SiC_osantowski.nk",
        "luxpop/SiC_palik.nk",
        "luxpop/SiC_windt.nk",
        "luxpop/SiC_yanagihara.nk",
        'scout/SiC (MIR).b',
        'scout/SiC (NIR-UV).b',
        #'scout_reverse_engineering/scout/SiC (simple infrared model).b'
        ]#}}}
## == SiO2 == #{{{
SiO2_files = [
        'scout/SiO2 (MIR-VUV).b',
        'scout/SiO2 [micron].b',
        'scout/SiO2 (fused).b',
        'sopra/SIO2.MAT',
        #'kitamura/SiO2/Koike1989.n',
        meep_materials.material_SiO2(where = None) ]
##SiO2_files += ['kitamura/SiO2/'+f for f in os.listdir('kitamura/SiO2/') if (f[-2:]=='.n' or f[-2:]=='.k')]
#}}}
## == Au == #{{{
Au_files = [
        meep_materials.material_Au(),
        'sopra/AU.MAT',
        'scout/Au (J, & C,, L, & H,).b',
        'scout/Au (JC).b',
        'scout/Au (MQ).b',
        'scout/Au [micron].b',
        'scout/Au.b',
        'scout/Au model.b',
        ]
#}}}
## == Si == #{{{
Si_files = [
        meep_materials.material_Si_NIR(),
        meep_materials.material_Si_MIR(),
        'other/Si_Dai2003.k',                  'other/Si_Dai2003.n',
        'sopra/SI100_2.MAT',
        'scout/Si (100).b',
        'scout/Si (Aspnes).b',
        'scout/Si (Vis-UV, Brendel model).b',
        'scout/Si (cryst,).b',
        'scout/Si (infrared).b',
        'scout/Si (crystalline, MIR-VUV).b',
        ]
#}}}
## == InP == #{{{
InP_files = [
        meep_materials.material_InP(),
        'sopra/INP.MAT',
        'scout/InP.b',
        'scout/InP (IR model).b',
        'scout/InP (Jellison).b',
        ]
#}}}
## == GaAs == #{{{
GaAs_files = [
        meep_materials.material_GaAs(),
        'sopra/GAAS.MAT',
        'sopra/GAAS031T.MAT',
        'scout/GaAs (100).b',
        'scout/GaAs.b',
        'scout/GaAs (31 deg C).b',
        ]
#}}}
## == InSb == #{{{
InSb_files = [
        'sopra/INSB.MAT',
        'scout/InSb (Jellison).b',
        'scout/InSb.b',
        ]
#}}}

## == STO == #{{{
STO_files = [
        'other/N_model_STO_300K.nk',
        'other/N_STO_300K.nk',
        'other/STO_Neville1972_300K.epsilon',
        #'other/STO_Neville1972_090K.epsilon',
        #'other/STO_Neville1972_004K.epsilon',
        #meep_materials.material_STO_THz(),
        meep_materials.material_STO(),
        ]
#}}}
## == Al2O3 == #{{{
Al2O3_files = [
        meep_materials.material_Sapphire(),
        meep_materials.material_Sapphire(ordinary=.66),
        meep_materials.material_Sapphire(ordinary=0),
        'other/N_Al2O3_c-cut.nk',
        'other/N_model_Al2O3_c-cut.nk',
        'scout/Al2O3.b',
        'scout/Al2O3 (Palik).b',
        'sopra/AL2O3.MAT',
        'sopra/AL2O3P.MAT',
        #'other/',
        ]
#}}}

## == Fe2O3 == #{{{
Fe2O3_files = [
        'scout/Fe2O3.b',
        ]
#}}}

## == Fe2O3 == #{{{
Fe3O4_files = [
        'scout/Fe3O4.b',
        ]
#}}}

## == Fe2O3 == #{{{
Cu_files = [
        'scout/Cu.b',
        ]
#}}}

## == Fe2O3 == #{{{
Ni_files = [
        'scout/Ni.b',
        ]
#}}}

## == Fe2O3 == #{{{
ZnO_files = [
        'scout/ZnO.b',
        ]
#}}}

## Dictionary of materials
materials = {
    'TiO2':     TiO2_files,
    'SiC':      SiC_files,
    'SiO2':     SiO2_files,
    'Au':       Au_files,
    'Si':       Si_files,
    'InP':      InP_files,
    'GaAs':     GaAs_files,
    'InSb':     InSb_files,
    'STO':      STO_files,
    'Al2O3':    Al2O3_files,
    'Fe2O3':    Fe2O3_files,
    'Fe3O4':    Fe3O4_files,
    'Ni':    Ni_files,
    'ZnO':    ZnO_files,
    'Cu':    Cu_files,
    }

##and not (len(sys.argv)>2 and (sys.argv[2] == "--big"))

if len(sys.argv)>1 and (sys.argv[1] in materials.keys()):
    data_sources    = materials[sys.argv[1]]
    output_file     = sys.argv[1]
else:
    print("Error: the first parameter may be one of the following materials:", materials.keys())
    exit(1)

if len(sys.argv)>2 and (sys.argv[2] == "--big"):
    output_file+="_big";  plt.figure(figsize=(16,16))
else:
    plt.figure(figsize=(8,8))

freq_range=(1e11, 1e16)
for number, data_source in enumerate(data_sources):
    try:
        if type(data_source)==str:
            plotlabel = data_source
            ## Select the file type to be loaded from
            n, eps, k = None, None, None
            print("Loading ", data_source)
            if data_source.endswith('.b'):        
                freq, eps = load_SCOUT_permittivity(data_source)
                n, k        = (eps**.5).real, (eps**.5).imag
                print(freq)
                print(eps)
                print(n)
                print(k)
            elif data_source.endswith('.MAT'):    
                freq, eps = load_SOPRA_permittivity(data_source)
                n, k        = (eps**.5).real, (eps**.5).imag
            elif data_source.endswith('.epsilon'):     
                freq, eps = np.loadtxt(data_source, usecols=[0,1], unpack=True, comments='#')
                n, k        = (eps**.5).real, (eps**.5).imag
                print(freq, n)
            elif data_source.endswith('.nk'):     
                freq, n, k, eps = load_n_k_eps(data_source)
            elif data_source.endswith('.n'):      
                lambda_micron, n = np.loadtxt(data_source, usecols=[0,1], unpack=True, comments='#')
                freq = 2.997e8 / (lambda_micron*1e-6)
                k, eps = None, None
            elif data_source.endswith('.k'):      
                lambda_micron, k = np.loadtxt(data_source, usecols=[0,1], unpack=True, comments='#')
                freq = 2.997e8 / (lambda_micron*1e-6)
                print(freq)
                n, eps = None, None
        ## Read the MEEP's material model
        if hasattr(data_source, 'pol'): 
            freq = 10**np.arange(np.log10(freq_range[0]), np.log10(freq_range[1]), .003)
            eps = analytic_eps(data_source, freq)
            n   = (analytic_eps(data_source, freq)**.5).real
            k   = (analytic_eps(data_source, freq)**.5).imag
            plotlabel = "model " + getattr(data_source, 'shortname', data_source.name)

        ## Plot the data
        if hasattr(data_source, 'pol'): 
            color = 'black'
        else:
            color = matplotlib.cm.hsv(float(number)/len(data_sources))

        if eps.all()!=None or (n.all()!=None and k.all()!=None):
            if eps.all()==None: eps=(n+1j*k)**.5
            plt.subplot(3,1,1)
            print("  Plotting epsilon for '%s' with %d data points" % (plotlabel, len(eps)))
            plt.plot(freq, eps+float(number)/100, color=color, marker='o', markersize=0, label=plotlabel)
            plt.plot(freq, eps.imag+float(number)/100, color=color, marker='s', markersize=0, ls='--')
            plt.xlim(freq_range); 
            plt.grid(True)
        if n.all()!=None:
            plt.subplot(3,1,2)
            print("  Plotting N for '%s' with %d data points" % (plotlabel, len(n)))
            plt.plot(freq, n, color=color, marker='o', markersize=0, label=plotlabel)
            plt.xlim(freq_range);
            plt.grid(True)
        if k.all()!=None:
            plt.subplot(3,1,3)
            print("  Plotting k for '%s' with %d data points" % (plotlabel, len(k)))
            plt.plot(freq, k, color=color, marker='o', markersize=0, label=plotlabel)
            plt.xlim(freq_range); 
            plt.grid(True)
    except:
        print("WARNING: data from file '%s' could not be plotted due to error" % (plotlabel))
        print(traceback.format_exc())

        print(sys.exc_info())
plt.subplot(3,1,1); plt.legend(loc=(0,0),prop={'size':6}); plt.ylabel("permittivity $\\epsilon_r'$"); plt.title(sys.argv[1])
plt.subplot(3,1,2); plt.legend(loc=(0,0),prop={'size':6}); plt.ylabel("index of refraction $n$");
plt.subplot(3,1,3); plt.legend(loc=(0,0),prop={'size':6}); plt.ylabel("index of absorption $k$");
plt.xlabel("Frequency $f$ [Hz]");
plt.show()
plt.savefig(output_file)

