# Case-cohort interval-censored data with time-dependent covariates
This SAS macro can be applied to interval-censored data with time-dependent covariates in a case-cohort study design. 

The work is built upon Sparling's SAS macro PARM_ICE.SAS (Version 3.2). By specifying a positive number (no greater than 1) after argument "prob=", user can indicate the sampling probability for a case-cohort design. "reps=" sets the number of resampling for bootstrapping. It is set to 10 by default. "sampleseed=" sets the initial seed for resampling. It is set to 0 by default. To initialize random number generation, a seed must be a positive integer. If you leave it as default or specify a negative seed, an initial seed is obtained from the the time of day from the computer's clock. 

For details of other parameters, we recommend users refer to Sparling's SAS macro documentation. 

## Example

Please check example_CCOH_PARM_ICE.SAS. 

## Reference: 

Gao X., Hudgens M., and Zou, F. (2021+), Case-cohort interval-censored data with time-dependentcovariates.

Sparling YH, Younes N, Lachin JM, and Bautista OM. Parametric survival models for interval-censored data with time-dependent covariates.Biostatistics;7, 599-614, 2006.
 
