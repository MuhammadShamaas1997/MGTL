Files:
	
plot_materials.py
	  - Reads and plots curves of permittivity, index of refraction or absorption.
	
plot_sSNOM_response.py
	  - using the permittivity function, plots the expected response of a s-SNOM microscope 
	
meep_materials.py
	  - library of materials for the MEEP simulation program 
	    [see http://fzu.cz/~dominecf/misc/meep/index.html]
	
scout
	sopra
	  - You may add the SOPRA or SCOUT database to get additional sources of data. 
	    I did not wish to redistribute the files for my page, although they may be 
	    freely downloaded from 
	    http://sspectra.com/sopra.html and http://www.wtheiss.com/download/scout.zip
	luxpop
	  - data downloaded from public web interface at http://luxpop.com/RefractiveIndexList_v2.html
	
kitamura
	  - data extracted from the http://www.seas.ucla.edu/~pilon/downloads.htm
For more information on obtaining material parameters, see
http://fzu.cz/~dominecf/misc/meep/index.html 


If available, both data from the aforementioned sources and from the Lorentzian model are used. 
This model is read from the meep_materials.py file, which is compatible with python-meep FDTD 
simulations.

See the source code of plot_materials.py for more information; 
run the little make.sh
file for examples of usage. Edit these files to add new materials or sources.

See also http://fzu.cz/~dominecf/misc/eps/index.html for the finished graphs of spectra


(c) Filip Dominec 2013
Licensed under BSD license terms
