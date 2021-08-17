*----------------------------------------------------------------------------------------------------------*;
* SAS Program    : example_CCOH_PARM_ICE.SAS                                                               *;
* Author         : Xiaoming Gao                                                        					   *;
* Abstract		 : This is an example of using the SAS macro CCOH_PARM_ICE.SAS                      	   *;
*                                                                                                          *;
* Reference      : 1) Gao X., Hudgens M., and Zou, F. (2021+), Case-cohort interval-censored data with     *;
*				   time-dependent covariates.  														   	   *;
*                  2) Sparling YH, Younes N, Lachin JM, Bautista OM. Parametric survival models for   	   *;
*                  interval censored data with time-dependent covariates.                                  *;
*                  Biostatistics (2006) 7(4):599-614.                                                      *;
*                                                                                                          *;
*                                                                                                          *;
* Date of Most Recent Revision : 16-Aug-2021                                                			   *;
*                                                                                                          *;
* SAS Version    : 9.4.                                         										   *;
*                                                                                                          *;
* This program is retrieved from: https://github.com/sophiegxm/CCOH_PARM_ICE				 			   *;
*----------------------------------------------------------------------------------------------------------*;
%let macpath=;/*Spedify location of the SAS macro*/
%let datapath=;/*Spedify location of the working data*/

/*----------------------------------------------------------------------------------------------------------*/
* The example dataset includes 500 subjects. Each subject has different number of visits and has two fixed  *;
* covariates x1, x2 and two time-dependent covariates u and y. Event=1 if a subject is interval-censored    *;
* between current and last visits. 																			*;
/*----------------------------------------------------------------------------------------------------------*/
LIBNAME loc "&datapath";
data analysis;
  set loc.exampledata;
run;
/*Windows*/
filename ccoh "&macpath\CCOH_PARM_ICE.SAS";
/*UNIX*/
/*filename ccoh "&macpath/CCOH_PARM_ICE.SAS";*/
%include ccoh;

Title "Simulated Data";
Title2 "Weibull Model Fit (Kappa=0)";
%CCOH_PARM_ICE(data=analysis,subject=id,t1=t1,t2=t2,fixed= X1 X2  ,
          tdctime=tau,TDC=  y u ,kappa=0 | 0, reps=2,prob=0.5,sampleseed = 20);
