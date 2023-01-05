# Dose_Interval
 
Code for the analysis presented in the scientific paper "The Impacts of SARS-CoV-2 Vaccine Dose Separation and Targeting on the COVID-19 epidemic in England", currently in press Nature Communications. An early pre-print is available at: https://doi.org/10.1101/2022.08.22.22278973.

This is a reduced version of the code used to generate the four model outputs (0=default, 1=3-week delay default vaccine efficacy, 2=3-week delay lower vaccine efficacy, 3=youngest first). 

The Matlab file Run_and_Plot.m runs the simulation code and generates graphs of the number of deaths and hospital admissions. For compactness this code only runs a single set of posterior values and for confidentiality reasons uses surrogate vaccination and population size data.

Run_and_Plot.m calls Generate_Output.m which in turn calls Simulate_One_Region.m; the underlying ODEs are within LeakyVacc_ODEs.m

The code can be made to run far more rapidly by using parfor loops and by converting LeakyVacc_ODEs.m and ODE_to_Observables.m into mex files.
