# Read Me for analysis code supporting 
"SOUTHERN OCEAN CO2 OUTGASSING AND NUTRIENT LOAD REDUCED BY A WELL-VENTILATED GLACIAL NORTH PACIFIC" 
Madison Shankle 
April 2025


## OVERVIEW: 
There are five python scripts performing the analysis supporting this study, described below. They should be contained within a directory called "code" and be run in succession. Before running, the "data" directory should exist at the same level as the "code" directory, containing the output from the four Earth system models (UVic, LOVECLIM, LOVECLIM-LGM, and GENIE). There should also exist at this level an empty directory  called "results", which is where generated figures will be saved. Finally, some scripts will print values to the console, such as regional averages and other values cited in the main text. Scripts 1 and 2 take ~3-5 minutes each to run, given the PCO2 and flux calculations they preform. All the other scripts should run quite quickly. 


## 01_NPacSObgc_calcPCO2_wPO4_andSigma.py
First analysis script. Calculates PCO2 from raw Earth system model (ESM) output. Prepares n=8 .nc files for use in subsequent scripts: "model_ctrl.nc" for control run or "model_pmoc.nc" for ventilated NPac, where "model" = GENIE, LC (LOVECLIM), LGM (LOVECLIM-LGM), or UV (UVic). Also generates the following figures: 4-paneled PCO2 anomaly figure for each ESM (Fig. 3) and same for PO43- anomalies (supp fig). Figure of Pacific depth-transect of PCO2 evolving over time (supp fig). 


## 02_NPacSObgc_moreFigs_wPO4_andSigma.py
Second analysis script. Generates additional figures related to PCO2 and calculates and prints regional and surface averages of PCO2 and PO43- anomalies cited in main text. Also calculates total air-sea CO2 flux from model output and decomposes it into anomalies due to various factors (temp and sal changes, wind changes, sea ice changes, etc.) Also generates the following figures: longitudinal transects of PCO2 and PO43- anomalies (Fig. 4), 5-panel figure of surface Southern Ocean anomalies (Fig. 5), air-sea CO2 flux anomaly decomposed into various factors (supp fig).  


## 03_NPacSObgc_StratWindsBio.py
Third analysis script. Performs calculations of changes in wind stress, stratification, and biological production. Calculates and prints regional and surface averages cited in main text. No figures generated.  


## 04_NPacSObgc_d14C_SF.py
Fourth analysis script. Calculates and generates figures of meridional stream function and deep-minus-surface radiocarbon age to illustrate and quantify induced ventilation (supp figs.)  


## 05_NPacSObgc_SurfFluxPCO2.py
Fifth analysis script. Performs additional analysis of surface pCO2 and air-sea CO2 flux. Takes zonal-integrals of both, generating integrated pCO2 and CO2 flux as a function of latitude. Prints some quantities of these. Generates figures: maps of absolute and anomaly PCO2 and flux (supp fig), 2-panel figure of zonally-integrated PCO2 and flux as a function of latitude (supp fig).  


## 06_NPacSObgc_biol.py
Sixth script. Generates a 6-paneled supplemental figure of global maps of surface phosphate concentration alongside net primary production rates, and their anomalies.  

