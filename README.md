# inhaled-RIF-popPK
sample-level ODE model for the comparison of an oral dose and inhaled dose of the anti-TB medication rifampin.
Based on the model presented in Ramachandran & Gadgil, 2023; DOI: https://doi.org/10.1002/psp4.13008
Lung absorption model based on Himstedt et al., 2022; DOI: https://doi.org/10.1093/jac/dkac240

Requires: Statistics and Machine Learning Toolbox, Parallel Computing Toolbox

## Main Scripts
### [`run_popPK_RIF`](./run_popPK_RIF)
Standard main script for the model. Calculates and plots concentration-time courses for a sample of patients for both an oral and inhaled dose of rifampin. Compares PK metrics (AUC and C_{max}), probability of target attainment (PTA), and cumulative fraction of response (CFR) for each method.

Key model parameters:
- n_pts: number of patients to simulate
- n_days_RIF: number of days to run the simulation for (set to > 3 for steady state)
- days_to_plot_RIF: number of days to include in the plots (backward from the last day in n_days)
- relevant_compts_RIF: compartments to plot and calculate PK metrics, PTA, and CFR for
- oral/lung_dose_RIF: oral dose amount in mg
- oral/lung_dose_freq_RIF: oral/lung dose frequency in doses/day

### optimize_regimen_RIF.m
Script to compare inhaled RIF dosing regimens, varying the dose amount and frequency for a combinatorial analysis. Compares all specified regimens of an inhaled dose to the standard oral dose of 600 mg 1x/day, then outputs a contour plot of the percent difference between the inhaled dosing regimen and the standard oral dose for AUC and C_{max}. Does not incorporate sample variation.

Key model parameters:
- n_days_RIF: number of days to run the simulation for
- n_days_RIF: number of days to run the simulation for (set to > 3 for steady state)
- relevant_compts_RIF: compartments to compare dosing regimens in
- lung_dose_min_RIF: minimum lung dose amount to simulate in mg
- lung_dose_max_RIF: maximum lung dose amount to simulate in mg
- lung_dose_inc_RIF: increment in which to vary the lung dose
- lung_dose_freq_min_RIF: minimum lung dose frequency to simulate in doses/day
- lung_dose_freq_max_RIF: maximum lung dose frequency to simulate in doses/day

## ODEs
The ODE systems for an oral dose model and an inhaled dose model are contained in oral_dose_ODEs and lung_dose_ODEs respectively, in the functions [`RIF_oral_ODEs.m`](./RIF_oral_ODEs.m) and [`RIF_lung_ODEs.m`](./RIF_lung_ODEs.m). These are based on the model in Ramachandran & Gadgil, with the absorption in RIF_lung_ODEs based on Himstedt et al. All compartment indices are listed at the top of each file, and all equation terms are labeled with the flow they represent.

## Methods
### getParamPDs.m
Returns cell arrays of probability distributions for each physiological parameter of the model based on coefficients of variation passed to it (one for all volume parameters, one for all flow parameters). All parameter distributions are normal and truncated at +/- 3 SDs.

Some model parameters are taken as fractions of body weight (for volume parameters) and cardiac output (for flow parameters). These are handled and returned separately from raw volume and flow parameters for renormalization later (so total volume/flow doesn't exceed 100% of body weight or cardiac output).

### loadPhysParams.m
Uses the probability distributions returned from get_param_PDs to sample the parameters of a specific patient. Renormalizes parameters taken as fractions of other model parameters. Returns two structs, one for physiological parameters (phys) and one for parition coefficients (pt - no variation incorporated).

### solveODEs.m
Takes the parameters for both RIF and the model and packages them to solve the model ODEs in RIF_oral/lung_ODEs for a specific patient. Returns matrices of concentration-time courses for each dosing method with each compartment as a column. These matrices are later packaged into cell arrays outside of the function.

### plotTimecourses.m
Plots the concentration-time courses returned by solveODEs. Returns both raw plots for all patients and plots of the 10th, 50th, and 90th percentiles for the patient sample. All figures are formatted for comparison as "oral plot | lung plot".

### trackPKMetrics.m
Calculates AUC and C_{max} from the concentration-time courses returned by solveODEs. Calculates the mean and SD for AUC and C_{max} for each dosing method for the patient sample. Conducts pairwise t-tests between dosing methods and returns the p-value and effect size. Writes output to popPK_analysis_day(n_days_RIF).xlsx.

### plotPTAs.m
Calculates and plots PTA and CFR for each nontoxic compartment based on given target values for AUC and C_{max} and a RIF-TB MIC distribution (modeled as a lognormal dist., based on data from the literature). Results are written to popPK_CFRs.xlsx.
