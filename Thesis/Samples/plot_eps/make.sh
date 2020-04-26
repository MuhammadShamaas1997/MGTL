#!/bin/bash
for material in SiC GaAs Si Au TiO2 SiO2 InP InSb Al2O3 STO
do 
    python ./plot_materials.py $material 
    python ./plot_materials.py $material --big
done

## The following is applicable to the s-SNOM microscopy simulation:
#for material in SiC GaAs Si Au TiO2 SiO2 InP      Al2O3   
#do 
#	python ./plot_sSNOM_response.py $material 
#	python ./plot_sSNOM_response.py $material --big
#done
