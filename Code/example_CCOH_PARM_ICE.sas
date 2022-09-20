*----------------------------------------------------------------------------------------------------------*;
* SAS Program    : example_CCOH_PARM_ICE.SAS                                                               *;
* Author         : Xiaoming Gao                                                        					   *;
* Abstract		 : This is an example of using the SAS macro ccoh_parm_ice.SAS                      	   *;
*                                                                                                          *;
* Reference      : 1) Gao X., Hudgens M., and Zou, F. (2022). Case-cohort studies with time-dependent      *;
*				   covariates and interval-censored outcome. Emerging Topics in Modeling Interval-Censored *;
*				   Survival Data (Editors J. Sun and D.G. Chen), ICSA Book Series in Statistics.  		   *;
*                  2) Sparling YH, Younes N, Lachin JM, Bautista OM. (2006). Parametric survival models for*;
*                  interval censored data with time-dependent covariates. Biostatistics 7(4):599-614.      *;
*                                                                                                          *;
*                                                                                                          *;
* Date of Most Recent Revision : 04-Mar-2022                                                			   *;
*                                                                                                          *;
* SAS Version    : 9.4.                                         										   *;
*                                                                                                          *;
* Github                                                                                                   *;
*----------------------------------------------------------------------------------------------------------*;
%let macpath=;/*Spedify location of the SAS macro*/
%let datapath=;/*Spedify location of the working data*/

LIBNAME loc "&datapath";
data ds;
  set loc.exampledata;
run;
/*Windows*/
filename ccoh "&macpath\ccoh_parm_ice.SAS";
/*Linux*/
/*filename ccoh "&macpath/ccoh_parm_ice.SAS";*/
%include ccoh;

Title "Simulated Data";
Title2 "Weibull Model Fit (Kappa=0)";
%CCOH_PARM_ICE(data=ds, subject=id,t1=t1,t2=t2,fixed= X1 X2  ,
	          tdctime=tau,TDC=y u,kappa=0 | 0,prob=0.2);
