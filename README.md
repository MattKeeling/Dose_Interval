# Dose_Interval
 
Code for the analysis presented in the scientific paper "The Impact of SARS-CoV-2 Vaccine Dose Separation and Dose Targeting on Hospital Admissions and Deaths from COVID-19 in England".

Preprint details: Matt J Keeling, Bridget Penman, Edward M Hill, Samuel Moore. (2022) The Impact of SARS-CoV-2 Vaccine Dose Separation and Dose Targeting on Hospital Admissions and Deaths from COVID-19 in England. *medRxiv*. doi: 10.1101/2022.08.22.22278973. URL: https://doi.org/10.1101/2022.08.22.22278973.

The Matlab file *Run_and_Plot.m* runs the simulation code and generates graphs of the impact of the modelled vaccine prioritisation scenarios. 
For compactness, this code only runs a single set of posterior values and uses surrogate vaccination and population size data.

*Vaccine_Compare.m* calls *Generate_Output.m* which in turn calls *Simulate_One_Region.m*; the underlying ODEs are within *LeakyVacc_ODEs.m*

The code can be made to run far more rapidly by using parfor loops and by converting *LeakyVacc_ODEs.m* and *ODE_to_Observables.m* into mex files.