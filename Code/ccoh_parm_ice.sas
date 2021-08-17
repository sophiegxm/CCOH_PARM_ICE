*----------------------------------------------------------------------------------------------------------*;
* SAS Program    : CCOH_PARM_ICE.SAS (Version 1.0)                                                         *;
* Function       : Fits an accelerated failure time model with time dependent covariates in a case-cohort  *;
*					design. 																			   *;
*                                                                                                          *;
* Reference      : 1) Gao X., Hudgens M., and Zou, F. (2021+), Case-cohort interval-censored data with     *;
*				   time-dependent covariates.  														   	   *;
*                  2) Sparling YH, Younes N, Lachin JM, Bautista OM. Parametric survival models for   	   *;
*                  interval censored data with time-dependent covariates.                                  *;
*                  Biostatistics (2006) 7(4):599-614.                                                      *;
*                                                                                                          *;
* Author         : Xiaoming Gao                                                        					   *;
* Date Completed : 08-Mar-2021 (Version 1.0)                                                               *;
*                                                                                                          *;
* This SAS macro is built on PARM_ICE.SAS (Version 3.2) written by Oliver M. Bautista                      *;
* Date of Most Recent Revision : 27-Jul-2021 (Version 1.0)                                                 *;
*                                                                                                          *;
* SAS Version    : This macro was developed using SAS version 9.4.                                         *;
* Printer Setup  : Font=SAS monospace 8pt ls=112 ps=60                                                     *;
*                  topmargin=0.5in bottommargin=0.5in leftmargin=0.5in rightmargin=0.5in                   *;
*                                                                                                          *;
* This macro is retrieved from: https://github.com/sophiegxm/CCOH_PARM_ICE				 				   *;
*                                                                                                          *;
*----------------------------------------------------------------------------------------------------------*;
%MACRO CCOH_PARM_ICE( DATA	 =                  /* SAS dataset containing the following variables:                        */
				,SUBJECT =                  /* Subject ID (must be numeric).                                          */
				,T1		 =                  /* Interval start time.                                                   */
				,T2		 =                  /* Interval end time.                                                     */
				,FIXED	 =                  /* List of fixed/baseline covariates.                                     */
				,TDCTIME =                  /* Time when the TDC value changed (Tau in the paper)                     */
				,TDC	 =                  /* List of time dependent covariates.                                     */
				,KAPPA	 = 0 | 0            /* Set the value to the left side of the delimiter equal to NULL/Blank to */
								            /* estimate an arbitrary distribution. Set the value to the right side    */
								            /* of the delimiter equal to the desired initial value to use in the      */
								            /* iterative estimation of KAPPA.                                         */
				,KAPPA_Ho = 0               /* Assign to this macro parameter a value corresponding to the null       */
								            /* hypothesis value of KAPPA, i.e., assumed null distribution. This       */
								            /* value will be used in the null model when KAPPA is to be estimated     */
								            /* by the macro. This value will also be used to do a LRT on the estimate */
								            /* of KAPPA when KAPPA is to estimated by the macro. The value assigned   */
								            /* to this parameter needs to be bounded by the limits specified in the   */
								            /* macro parameter KAPPA_LIM. The value assigned to this parameter is     */
								            /* irrelevant when KAPPA is held fixed through the macro parameter KAPPA. */
				,KAPPA_LIM = -3 | 3         /* Set the value to the left side of the delimiter equal to the lower     */
								            /* limit of KAPPA, the right side of the delimiter equal to the upper     */
								            /* limit of KAPPA that will be used in the iterative estimation of KAPPA  */
								            /* when KAPPA is not held fixed.                                          */
				,EPSILON   = 0.01           /* Specify a number between 0 and 1.0                                     */
				,SCALE     =   | 1          /* Set the value to the left side of the delimiter equal to NULL/Blank to */
								            /* estimate an arbitrary value of the SCALE parameter. Set the value to   */
								            /* the right side of the delimiter equal to the desired initial value to  */
								            /*  use in the iterative estimation of the SCALE parameter.               */
				,SCALE_LB  = 0.00000000001  /* Lower bound of the SCALE parameter.                                    */
				,DETAILS   =                /* Valid values are NON, SHORT, or LONG.                                  */
				,SAVELIB   = WORK           /* Name of library where COVPARMS, PARMEST, and LRTMODEL will be saved.   */ 
				,DELIMITER = |              /* Delimiter character used in specifying fixed and initial KAPPA values. */ 
				,prob =               		/* Probability of inclusion into subcohort.                               */ 
				,reps = 10              	/* Number of resampling for bootstrapping.                                */ 
				,sampleseed = 0             /* Seed for resampling                                                    */ 
		);
 *---------------------------------------------------------------------------------------------------------*;
 * Notes on T1 and T2 (M=missing, NM=non-missing)                                                          *;
 * ---------------------------------------------                                                           *;
 *   1. T1=T2(NM)     : Exact event time                                                                   *;
 *   2. T1(M) ,T2(NM) : Left censored observation                                                          *;
 *   3. T1(NM),T2(M)  : Right censored observation                                                         *;
 *   4. T1(NM),T2(NM) : Interval censored observation                                                      *;
 *   5. T1(M) ,T2(M)  : invalid observation that the program will delete.                                  *;
 *                                                                                                         *;
 * Notes on KAPPA                                                                                          *;
 * --------------                                                                                          *;
 *   KAPPA = blank (arbitrary distn, i.e. unrestricted Kappa)                                              *;
 *         = 0 (Weibull)                                                                                   *;
 *         = 1 (Log-Logistic)                                                                              *;
 *                                                                                                         *;
 * Notes on Epsilon                                                                                        *;
 * ----------------                                                                                        *;
 *  EPSILON is a number between 0 and 1.0 such that when (KAPPA-1)< EPSILON, the form of the hazard and    *;
 *  gradient functions corresponding to the case when KAPPA=1 will be used.                                *;
 *                                                                                                         *;
 * Notes on DETAILS                                                                                        *;
 * ----------------                                                                                        *;
 * DETAILS = None  (supresses all printouts generated by the IML optimization routines).                   *;                                           *;
 *         = Short (prints the IML summary of optimization start, termination, and iteration history)      *;
 *         = Long  (prints initial and final parameters in addition to all printed by the Short option     *;
 *---------------------------------------------------------------------------------------------------------*;


   *-------------------------------------------------------*;
   * Section 1:   Extract Names of Fixed and TD Covariates *;
   *-------------------------------------------------------*;
   %Let N=0;
   %Do %While(%Scan(&FIXED,&N+1)^=);
     %Let N=%Eval(&N+1);
     %Let Var&N=%Scan(&FIXED,&N);
   %End;
   %Let M=0;
   %Do %While(%Scan(&TDC,&M+1)^=);
     %Let M=%Eval(&M+1);
     %Let L=%Eval(&M+&N);
     %Let Var&L=%Scan(&TDC,&M);
   %End;
   *-------------------------------------------------------------------------------------------------------------*;
   * Section 2:   Pre-processing of the input dataset.                                                           *;
   *-------------------------------------------------------------------------------------------------------------*;
   * A temporary analysis SAS dataset will be created containing the data that will be analyzed.                 *;
   * A.  Subjects with only 1 record will be flagged.                                                            *;
   *      a.1) In version 3.0 of the macro, single-record subjects will be excluded from the analysis dataset.   *;
   * B.  Records that satisfy the following conditions will be excluded from the analysis dataset:               *;
   *      b.1) Records with missing time interval start and end times;                                           *;
   *      b.2) Records with time interval start or end times that are negative;                                  *;
   *      b.3) Records with time interval start time that is greater than the end time.                          *;
   * C.  Check will be conducted to ensure that subjects have the minimum number of records.                     *;
   *     If the minimum number of records is not met, additional records will be created as follows:             *;
   *      c.1) Subjects with exact event times, left censored, or right censored should have at least 2 records. *;
   *      c.1.1) The first record will normally have TDCTIME=0 representing the baseline value of the TDC.       *;
   *      c.1.2) There should be a record where TDCTIME = max(T1,T2).                                            *;
   *      c.1.3) If there is no record corresponding to TDCTIME = 0, a record will be created. The value of the  *;
   *             TDC will be imputed by first observation carried backwards (i.e., TDC WILL BE IMPUTED BY FOCB). *;
   *      c.1.4) If there is no update of the TDC that occured at time = max(T1,T2), a record will be created    *;
   *             where TDCTIME = max(T1,T2) and the last observed value of TDC will be assigned to this created  *;
   *             record (i.e., TDC WILL BE IMPUTED BY LOCF method).                                              *;
   *      c.2) Subjects who are interval censored should have at least 3 records.                                *;
   *      c.2.1) The first record will normally have TDCTIME=0 representing the baseline value of the TDC.       *;
   *      c.2.2) There should be another record where TDCTIME = T1.                                              *;
   *      c.2.3) There should be another record where TDCTIME = T2.                                              *;
   *      c.2.4) If there is no record corresponding to TDCTIME = 0, a record will be created. The value of the  *;
   *             TDC will be imputed by first observation carried backwards (i.e., TDC WILL BE IMPUTED BY FOCB). *;
   *      c.2.5) If there is no update of the TDC that occured at time T1 and/or T2 record(s) will be created    *;
   *             such that there will be record(s) where TDCTIME = T1 and TDCTIME = T2. The value of TDC in      *;
   *             the record(s) that will be created will be IMPUTED by LOCF method.                              *;
   *-------------------------------------------------------------------------------------------------------------*;
   	 proc sql;
	 	create table InputDat as 
		select &subject., &tdctime., %do i = 1 %to &L.; &&Var&i , %end; &t1., &t2.,
			   case when &T1  = &T2				      then 0	/* Exact event time  */
         		    when (&T1 <= .Z) and (&T2 > 0)	  then 1	/* Left  censored    */
				    when (&T2 <= .Z) and (&T1 > 0)	  then 2	/* Right censored    */
				    when (0  <= &T1) and (&T1 < &T2)  then 3	/* Interval censored */
				    end 
				    as _censor_ LABEL='0=Exact, 1=Left, 2=Right, 3=Interval',
				count(&tdctime.) as nobs LABEL='Subject level - number of records'
		from &DATA.
		group by &subject.
		order by &subject., &tdctime.
		;
	quit;
	data ReducedData;
        set InputDat(where=(NOBS>1));                                    /* Exclude subjects based on a.1) */
		if &T1 <= .Z and &T2 <= .Z then delete;                          /* Exclude subjects based on b.1) */
        else if ((.Z < &T1 < 0) | (.Z < &T2 < 0)) then delete;           /* Exclude subjects based on b.2) */
        else if ((&T1 > .Z) and (&T2 > .Z) and (&T1 > &T2)) then delete; /* Exclude subjects based on b.3) */
	run;
	*---------------------------------------------------------*;
	* Add records as needed -- see item C of Section 2 header *;
	*---------------------------------------------------------*;
	*--------------------------------------*;
	*-- Check for presence of TDC time=0 --*;
	*--------------------------------------*;
	data Check_T0; 
		set ReducedData;
		by &subject.;
		output;
		if first.&subject. then do;
			if &tdctime. > 0 then do;
				&tdctime. = 0;
				output;
			end;
		end;
	run;
	proc sort data=Check_T0;
		by &subject. &tdctime.;
	run;
	*--------------------------------------------------------------*;
	*-- Check for presence of TDC time = event or censoring time --*;
	*--------------------------------------------------------------*;
	data censored_int censored_oth;
		set Check_T0;
		if _censor_ in (0,1,2) then output censored_oth;
		if _censor_ in (3)     then output censored_int;
	run;
	*--------------------------------------------------------------------------*;
	*-- Checks related to left censored, right censored, or exact event time --*;
	*--------------------------------------------------------------------------*;
	data censored_oth;
		set censored_oth;
		by &subject.;
		if last.&subject. then do;
			output;
			if &tdctime. < max(&t1.,&t2.) then do;
				&tdctime. = max(&t1.,&t2.); output;
			end;
		end;
		else output;
	run;
	*---------------------------------------------------*;
	*-- Checks related to interval censored endpoints --*;
	*---------------------------------------------------*;
		data dummy1;
			set censored_int;
			by &subject. ;
			if first.&subject. then do;
				&tdctime. = &t1. ; output;
				&tdctime. = &t2. ; output;
			end;
			drop %do i = %eval(&N.+1) %to &L.; &&Var&i %end;;
		run;
		data dummy2;
			merge dummy1 censored_int ;
			by &subject. &tdctime. ;
		run;
		data censored_int(drop = %do i = %eval(&N.+1) %to &L.; var_&i %end;);
			set dummy2;
			by &subject. ;
			retain %do i = %eval(&N.+1) %to &L.; var_&i %end; ;
			if first.&subject. then do;
				%do i = %eval(&N.+1) %to &L.;
					var_&i = &&Var&i ;
				%end;
			end;
			else do;
				if ((&tdctime.=&t1.) or (&tdctime.=&t2.)) then do;
					%do i = %eval(&N.+1) %to &L.;
						if &&Var&i = . then &&Var&i = var_&i. ;
						else var_&i = &&Var&i ;
					%end;
				end;
				else do;
					%do i = %eval(&N.+1) %to &L.;
						var_&i = &&Var&i ;
					%end;
				end;
			end;
		run;
	*-------------------------------------------*;
	*-- Combine censored_int and censored_oth --*;
	*-------------------------------------------*;
	data Recs_Added;
		set censored_oth censored_int;
	run;
     *---------------------------------------------------*;
     * Get arguments passed to the macro parameter KAPPA *;
     *---------------------------------------------------*;
		%if %index(&kappa.,&delimiter.)=1 %then %do;
			%LET FIXED_KAPPA = %scan(&kappa.,-2,"&delimiter.");
			%LET INIT_KAPPA  = %scan(&kappa., 1,"&delimiter.");
		%end;
		%if %index(&kappa.,&delimiter.)>1 %then %do;
			%LET FIXED_KAPPA = %scan(&kappa., 1,"&delimiter.");
			%LET INIT_KAPPA  = %scan(&kappa., 2,"&delimiter.");
		%end;
		%if &INIT_KAPPA =  %then %LET INIT_KAPPA = 0;
     *-------------------------------------------------------*;
     * Get arguments passed to the macro parameter KAPPA_LIM *;
     *-------------------------------------------------------*;
		%LET KAPPA_LB = %scan(&kappa_lim.,1,"&delimiter.");
		%LET KAPPA_UB = %scan(&kappa_lim.,2,"&delimiter.");
     *---------------------------------------------------*;
     * Get arguments passed to the macro parameter SCALE *;
     *---------------------------------------------------*;
		%if %index(&scale.,&delimiter.)=1 %then %do;
			%LET FIXED_SCALE = %scan(&scale.,-2,"&delimiter.");
			%LET INIT_SCALE  = %scan(&scale., 1,"&delimiter.");
		%end;
		%if %index(&scale.,&delimiter.)>1 %then %do;
			%LET FIXED_SCALE = %scan(&scale., 1,"&delimiter.");
			%LET INIT_SCALE  = %scan(&scale., 2,"&delimiter.");
		%end;
		%if &INIT_SCALE =  %then %LET INIT_SCALE = 1;

   PROC IML;   ** WORKSIZE=250000;
     *----------------------------------------------*;
     * Define Log-likelihood and Gradient functions *;
     *----------------------------------------------*;
     *-------------------------------------------------------------*;
     * ***** Log-likelihood function *****                         *;
     * The argument W to the function is the parameter vector at   *;
     * which the likelihood is calculated.                         *;
     * The argument is assumed partitioned as follows:             *;
     * w(|1|)=alpha                                                *;
     * w(|2|)=theta                                                *;
     * w(|3|)=kappa                                                *;
     * w(|4...p+3|)=gamma                                          *;
     * w(|p+4...p+q+3|)=eta                                        *;
     *-------------------------------------------------------------*;
     START LOGLIK(W) GLOBAL(N,P,Q,Z,K,T,TDC,CENSOR,START,STOP,EPSILON);
        Alpha = w(|1|);                 /* extract a few things to */
        Theta = w(|2|);                 /* ... simplify the code   */
        Kappa = w(|3|);
        Gamma = w(|4:(p+3)|);
        Eta   = w(|(p+4):(p+q+3)|);
        B0    = Theta + (Z*Gamma);      /* Fixed Part of Beta      */
     /*------------------------------------------------------------*/
     /*  Initialize Loglik                                         */
     /*------------------------------------------------------------*/
         LogLik=0;
           /*------------------------------------------------------*/
           /*  Compute CapL1, CapL2 and contribution to Log-      -*;
           /*  likelihood, by Subject ID.                         -*;
           /*------------------------------------------------------*/
            do i=1 to N;
                 CENS=CENSOR(|i|); T1=T(|i,1|); T2=T(|i,2|); Ki=k(|i|);
               *---------------------------------------------*;
               * Get TAU and TDC(i.e. Y) of subject i        *;
               *---------------------------------------------*;
                     TAU=TDC(|start(|i|):stop(|i|),2|);
                       Y=TDC(|start(|i|):stop(|i|),3:(2+Q)|);
               *---------------------------------------------*;
               * Compute subject contribution to cum hazard  *;
               *---------------------------------------------*;
				Log_Bi	  = B0(|i|)+(Y*Eta);
                B_i		  = exp(Log_Bi);
				DummyTAU  = TAU; DummyTau(|LOC(TAU=0)|)=.;
				LogTAU	  = log(DummyTAU); LogTAU(|LOC(TAU=0)|)=0;
				TAU_A	  = exp(LogTAU#Alpha); TAU_A(|LOC(TAU=0)|)=0;
				BTAU_R	  = (B_i#TAU_A)(|1:(ki-1)|);
				BTAU_L	  = B_i(|1:(ki-1)|)#TAU_A(|2:ki|);
				ONE_KAPPA = 1-Kappa;
                if Ki>1 then do;
                 if abs(ONE_Kappa) >  Epsilon then do;
						num_R = exp(log(1 + BTAU_R)#ONE_Kappa);
						num_L = exp(log(1 + BTAU_L)#ONE_Kappa);
						Right = num_R/ONE_Kappa;
						Left  = num_L/ONE_Kappa;
                 end;
				 if abs(ONE_Kappa) <=  Epsilon then do;
                    Right=log(1+BTAU_R);
                     Left=log(1+BTAU_L);
				 end;
                   Pi = ( 1 - (1/(1+BTAU_L)) )(|ki-1|);
                   	CapLUpdt=(Left-Right);
                   	TAU_High=TAU(|2:ki|);
                   CapL1=sum( CapLUpdt#(TAU_High<=T1) );
                   CapL2=sum( CapLUpdt#(TAU_High<=T2) );
                end;
                else do;
                   CapL1=0; CapL2=0; Pi=0;
                end;  *** End ki > 1 loop **;
                /*-----------------------------------------------*/
                /*  Calculate components of Log-likelihood       */
                /*-----------------------------------------------*/
                if Cens=0 then do;
					Lambda = alpha*Pi*((1-Pi)**(Kappa-1))/T1;
					if Lambda <= 0
						then Li = -capL1;
						else Li = Log(Lambda)-capL1;
                end;
                if Cens=1 then do;
					S2 = exp(-capL2);
					if 1- S2 <= 0
						then Li = 0;
						else Li = Log(1-S2);
                end;
                if Cens=2 then do;
					Li = -capL1;
                end;
                if Cens=3 then do;
					S1 = exp(-capL1);
					S2 = exp(-capL2);
					if S1-S2 <= 0
						then Li = 0;
						else Li = Log(S1-S2);
                end;
				/*Incorporate case-cohort weight*/
				if Cens=2 then weight = 1/&prob;
				else weight=1;
				Li = Li*weight;
                LogLik = Loglik + Li;**Sum all contributions to Loglik**;
            end; *** End i=1 to N do loop ***;
            return(LogLik);
     FINISH LOGLIK;
     *------------------------------------------------------------*;
     * ***** Gradient function *****                              *;
     * This function works exactly like the log-likelihood        *;
     * function, except more things are calculated.               *;
     * The argument W to the function is similar to the W in the  *;
     * Log-likelihood function.                                   *;
     * The argument is assumed partitioned as follows:            *;
     * w(|1|)=alpha                                               *;
     * w(|2|)=theta                                               *;
     * w(|3|)=kappa                                               *;
     * w(|4...p+3|)=gamma                                         *;
     * w(|p+4...p+q+3|)=eta                                       *;
     *------------------------------------------------------------*;
     START GRADIENT(W) GLOBAL(N,P,Q,R,Z,K,T,TDC,CENSOR,START,STOP,EPSILON);
        Grad  = J(r,1,0);              /* Initialize Gradient to 0 */
        Alpha = w(|1|);                /* Extract a few things to  */
        Theta = w(|2|);                /* simplify the code.       */
        Kappa = w(|3|);
        Gamma = w(|4:(p+3)|);
        Eta   = w(|(p+4):(p+q+3)|);
        B0    = Theta + (Z*Gamma);     /* Fixed Part of Beta       */
     /*------------------------------------------------------------*/
     /*  Initialize Gradient vector                                */
     /*------------------------------------------------------------*/
       Gradient=J(1,r,0); Grad1=Gradient; Grad2=Grad1; Grad_i=Grad2;
           /*------------------------------------------------------*/
           /*  Compute CapL1, CapL2, Grad1, & Grad2 by Subject ID  */
           /*------------------------------------------------------*/
			DO i=1 to N;
                 CENS=CENSOR(|i|); T1=T(|i,1|); T2=T(|i,2|); Ki=K(|i|);
                *---------------------------------------------*;
                * Get TAU and TDC(i.e. Y) of subject i        *;
                *---------------------------------------------*;
                     TAU=TDC(|start(|i|):stop(|i|),2|);
                       Y=TDC(|start(|i|):stop(|i|),3:(2+Q)|);
               *---------------------------------------------*;
               * Compute subject contribution to cum hazard  *;
               *---------------------------------------------*;
				Log_Bi	  = B0(|i|)+(Y*Eta);
                B_i		  = exp(Log_Bi);
				DummyTAU  = TAU; DummyTau(|LOC(TAU=0)|)=.;
				LogTAU	  = log(DummyTAU); LogTAU(|LOC(TAU=0)|)=0;
				TAU_A	  = exp(LogTAU#Alpha); TAU_A(|LOC(TAU=0)|)=0;
				BTAU_R	  = (B_i#TAU_A)(|1:(ki-1)|);
				BTAU_L	  = B_i(|1:(ki-1)|)#TAU_A(|2:ki|);
				ONE_KAPPA = 1-Kappa;
                if Ki>1 then do;
                 if abs(ONE_Kappa) >  Epsilon then do;
						num_R = exp(log(1 + BTAU_R)#ONE_Kappa);
						num_L = exp(log(1 + BTAU_L)#ONE_Kappa);
						Right = num_R/ONE_Kappa;
						Left  = num_L/ONE_Kappa;
                 end;
				 if abs(ONE_Kappa) <=  Epsilon then do;
                    Right=log(1+BTAU_R);
                     Left=log(1+BTAU_L);
				 end;
				 CapLUpdt=(Left-Right);
				 TAU_High=TAU(|2:ki|);
                *----------------*;
                * Tau[j+1] <= T1 *;
                *----------------*;
                 CapL1 = sum( CapLUpdt#(TAU_High<=T1) );
                *----------------*;
                * Tau[j+1] <= T2 *;
                *----------------*;
                 CapL2 = sum( CapLUpdt#(TAU_High<=T2) );
                *---------------------------------------------*;
                * Pi1={Pijp1}={ 1-1/(1+Bj*Tau[j+1]**alpha)}   *;
                * Pi0={Pij}  ={ 1-1/(1+Bj*Tau[j]**alpha)}     *;
                *---------------------------------------------*;
                  Pi1 = ( 1 - (1/(1+BTAU_L)) );
                  Pi0 = ( 1 - (1/(1+BTAU_R)) );
                *---------------------------------------------*;
                * Last Pi1 value                              *;
                *---------------------------------------------*;
                     Pi=Pi1(|ki-1|);
                *---------------------------------------------*;
                *            Derivative wrt Alpha             *;
                *---------------------------------------------*;
                   LTAU     = LogTau(|1:ki-1|);
                   LTAU1    = LogTau(|2:ki|);
 
				   if abs(1-Kappa) >  Epsilon then do;
	                   dA_Right = Right#Pi0#LTAU#(1-Kappa);
    	               dA_Left  = Left#Pi1#LTAU1#(1-Kappa);
				   end;
 
				   if abs(1-Kappa) <=  Epsilon then do;
	                   dA_Right = Pi0#LTAU;
    	               dA_Left  = Pi1#LTAU1;	               
				   end;
 
                   dA_Updt  = dA_Left - dA_Right;
                *---------------------------------------------*;
                *    Derivative wrt Beta-related quantities   *;
                *    (i.e. Theta, Gamma, Eta)                 *;
                *---------------------------------------------*;
				   if abs(1-Kappa) >  Epsilon then do;
						dB_Right = Right#Pi0#(1-Kappa);
						dB_Left  = Left#Pi1#(1-Kappa);
				   end;
				   if abs(1-Kappa) <=  Epsilon then do;
						dB_Right = Pi0;
						dB_Left  = Pi1;
				   end;
 
                   dT_Updt  = dB_Left - dB_Right;     *-- Theta --*;
                   dummy    = dt_Updt;
                   if max((dummy=.))=1 then do;
                      dt_Updt(|LOC((dummy=.))|)=0;
                   end;
                   dG_Updt  = (dT_Updt)*Z(|i,|) ;     *-- Gamma --*;
                   dE_Updt  = (dT_Updt)#Y(|1:ki-1,|); *-- Eta   --*;
                *---------------------------------------------*;
                *    Derivative wrt Kappa                     *;
                *---------------------------------------------*;
				   if abs(1-Kappa) >  Epsilon then do;
						Pi0m=(1-Pi0);
						if min(Pi0m)<=0 then do;
							DumPi0m=Pi0m; DumPi0m(|LOC(Pi0m<=0)|)=.;
							LogPi0m=Log(DumPi0m); LogPi0m(|LOC(Pi0m<=0)|)=0;
						end;
						else LogPi0m=Log(Pi0m);
						dK_Right = Right#(LogPi0m+(1/(1-Kappa)));
 
						Pi1m=(1-Pi1);
						if min(Pi1m)<=0 then do;
							DumPi1m=Pi1m; DumPi1m(|LOC(Pi1m<=0)|)=.;
							LogPi1m=Log(DumPi1m); LogPi1m(|LOC(Pi1m<=0)|)=0;
						end;
						else LogPi1m=Log(Pi1m);
						dK_Left  = Left#(LogPi1m+(1/(1-Kappa)));
						dK_Updt  = dK_Left - dK_Right;
					end;
 
					if abs(1-Kappa) <=  Epsilon then do;
						dK_Updt = 0;
					end;
                ******************;
                * Tau[j+1] <= T1 *;
                ******************;
                *--- Alpha ---*;
                   Grad1(|1|) = sum( dA_Updt#(TAU_High<=T1) );
                *--- Theta ---*;
                   Grad1(|2|) = sum( dT_Updt#(TAU_High<=T1) );
                *--- Kappa ---*;
                   Grad1(|3|) = sum( dK_Updt#(TAU_High<=T1) );
                *--- Gamma ---*;
                   Grad1(|4:p+3|)=(dG_Updt#(TAU_High<=T1))(|+,|);
                *--- Eta   ---*;
                   Grad1(|p+4:r|)=(dE_Updt#(TAU_High<=T1))(|+,|);
                ******************;
                * Tau[j+1] <= T2 *;
                ******************;
                *--- Alpha ---*;
                   Grad2(|1|) = sum( dA_Updt#(TAU_High<=T2) );
                *--- Theta ---*;
                   Grad2(|2|) = sum( dT_Updt#(TAU_High<=T2) );
                *--- Kappa ---*;
                   Grad2(|3|) = sum( dK_Updt#(TAU_High<=T2) );
                *--- Gamma ---*;
                   Grad2(|4:p+3|)=(dG_Updt#(TAU_High<=T2))(|+,|);
                *--- Eta   ---*;
                   Grad2(|p+4:r|)=(dE_Updt#(TAU_High<=T2))(|+,|);
			END;
			ELSE DO;
                   CapL1=0; CapL2=0; Pi=0;
                   Grad1(|1:r|)=0; Grad2(|1:r|)=0;
			END; ** End of k(|i|) > 1 do loop **;
               /*------------------------------------------------*/
               /*  Calculate components of Gradient vector       */
               /*------------------------------------------------*/
                if Cens=0 then do;
                  *--- Alpha ---*;
					DummyT1 = T1; DummyT1(|LOC(T1=0)|)=.;
					LogT1   = log(DummyT1); LogT1(|LOC(T1=0)|)=0;
                    Grad_i(|1|)=(1/alpha)+(1-kappa*Pi)*LogT1;
                  *--- Theta ---*;
                    Grad_i(|2|)=(1-kappa*Pi);
                  *--- Kappa ---*;
                    If Pi<1 then Grad_i(|3|)= Log(1-Pi);
                    else Grad_i(|3|)= 0;
                  *--- Gamma ---*;
                    Grad_i(|4:p+3|)=Grad_i(|2|)#Z(|i,|);
                  *--- Eta   ---*;
                    Grad_i(|p+4:r|)=(1-kappa*Pi)#Y(|ki-1,|);
                  *--- Gradient --*;
                    Gradient = Gradient + ( Grad_i - Grad1 );
                end;
                if Cens=1 then do;
					S2 = exp(-capL2);
					if 1-S2 ^= 0 
						then Gradient = Gradient + ((S2/(1-S2))#Grad2);
						else Gradient = Gradient;
                end;
                if Cens=2 then do;
                    Gradient = Gradient - Grad1/&prob; /*weight applied*/
                end;
                if Cens=3 then do;
					S1 = exp(-capL1);
					S2 = exp(-capL2);
					if S1-S2 ^= 0 
						then Gradient = Gradient+((S2#Grad2)-(S1#Grad1))/(S1-S2);
						else Gradient = Gradient;
                end;
               /*--------------------------------------------------*/
               /* End computation Gradient Computation.            */
               /*--------------------------------------------------*/
            end; ** End of i=1 to N do loop **;
        return(Gradient);
     FINISH GRADIENT;

	 /* store modules */
	 store module=LOGLIK;
	 store module=GRADIENT;
	quit;
	*----------------------------------------------------------*;
	* Set up Bootstrapping Datasets*;
	*----------------------------------------------------------*;

	proc sort data=Recs_Added;by &subject. ; run;
	data Recs_Added;
	  set Recs_Added;
	  by &subject.;
	  retain subject_N 0;
	  if first.&subject. then subject_N +1;
	  if last.&subject. then call symputx('nsub',subject_N);
	run;

	proc sql;
		create table ds_id as 
			select distinct &subject.,subject_N
			from Recs_Added;
	quit;

	/*------------------*/
	/* Resample subjects*/
	/*------------------*/
	PROC SURVEYSELECT DATA=ds_id OUT=rsample METHOD=URS
	  samprate = 1 reps=&reps SEED=&sampleseed noprint; 
	RUN;

	proc sort data=rsample;by &subject. subject_N Replicate; run;
	proc transpose data=rsample out=widersample prefix=rep;
		id replicate;
		by &subject. subject_N;
		var NumberHits;
	run;

	proc sort data=Recs_Added;by &subject. subject_N; run;
	data interval_bs ;
		merge Recs_Added(in=s) widersample(in=r);
		by &subject. subject_N;

		if r then output interval_bs;
	run;

	%DO nreps=1 %TO &reps;
	  data interval_rp;
		set interval_bs(where=(rep&nreps>.));;
		by &subject. &tdctime.;
		resample=&nreps;
		do i=1 to rep&nreps;
			&subject._r = 10**ceil(log10(&nsub))*i+subject_N; /*NEED TO BE NUMERIC*/
			output;
		end;
		drop i rep:;
	  run;

	  proc append base=interval_hits data=interval_rp force;
	  run;
	%END;

	data interval_hits2;
	  set interval_hits Recs_Added(in=r);
	  if r then do; resample=0; &subject._r =subject_N;end;
	run;
	/*-------------------*/
	/* End of Resampling*/
	/*-------------------*/
	*----------------------------------------------------------*;
	* End of Bootstrapping Datasets Set up*;
	*----------------------------------------------------------*;


	%DO nreps=0 %TO &reps;

		proc sql;
		 	create table UseData as 
			select &subject._r, &tdctime., %do i = 1 %to &L.; &&Var&i , %end; &t1., &t2., _censor_,
					count(&tdctime.) as K LABEL='Subject level - number of records'
			from interval_hits2 where resample=&nreps.
			group by &subject._r
			order by &subject._r, &tdctime.
			;
	     *-------------------------------------------------------------*;
	     * Create a dataset containing TDC and time of change of TDC   *;
	     *-------------------------------------------------------------*;
			create table TDCVARS as 
			select &subject._r ,&TDCTIME. %do i = %eval(&N.+1) %to &L.; ,&&Var&i %end;
			from UseData
			;
	     *-----------------------------------------------------------------------------------------------*;
	     * Create a dataset containing fixed covariates and event/censoring times - 1 record per subject *;
	     *-----------------------------------------------------------------------------------------------*;
			create table SUBJECT as 
			select DISTINCT &subject._r ,&T1. ,&T2. , _CENSOR_ %do i = 1 %to &N.; ,&&Var&i %end; ,K
			from UseData
			;
		 quit;
	     data SUBJECT;
	        set SUBJECT;
	        retain _START_ _STOP_;
	        if _N_=1 then do;
				_START_=1; _STOP_=K;
			end;
			else do;
				_START_ = _STOP_ + 1;
				_STOP_  = _STOP_ + K;
			end;
	        Label _START_="Start location of subject TDC"
	               _STOP_="Stop location of subject TDC";
	     run;
	     *---------------------------------------------------*;
	     * Add weight variable to each subject *;
	     *---------------------------------------------------*;
		 data UseData;
			set UseData;
			if _censor_=2 then weight = 1/&prob;
				else weight=1;
		 run;


	     *---------------------------------------------------------------------------------------------------*;
	     * Get initial value of the parameter vector Omega from a Weibull model, ignoring the time-dependent *;
	     * nature of the time-dependent covariates.                                                           *;
	     *---------------------------------------------------------------------------------------------------*;
		 proc sort data=usedata;by &subject._r  &tdctime;run;

		 data newuse;
			set usedata;
			by &subject._r  &tdctime;
			if last.&subject._r ;
		 run;

/*options errors=0;*/
	     proc lifereg data=newuse OUTEST=INITPARM noprint;
		 	weight weight;
	         model (&T1.,&T2.)=&FIXED &TDC / Distribution=weibull ;
	     run;
/*options errors=20;*/
		/* check if there is any observation*/
		data _NULL_;
		if 0 then set initparm nobs=n;
		call symputx('nrows',n);
		stop;
		run;
		/* Use exponential distribution if the model does not converge using Weibull*/

		 %put &nrows;
		 %if &nrows=0 %then %do; 
		     proc lifereg data=newuse OUTEST=INITPARM noprint;
			 	weight weight;* revision 07/18/2020;
		         model (&T1.,&T2.)=&FIXED &TDC / Distribution=EXP ;
		     run;
			 %if %symexist(lifereg_err)=0 %then %let lifereg_err=0;
			 	%else %let lifereg_err=%eval(&lifereg_err+1);

			%put &lifereg_err;
		 %end;	 
	   *--------------------*;
	   * Estimate TDC Model *;
	   *--------------------*;
	   PROC IML;   ** WORKSIZE=250000;
		load module=LOGLIK;
		load module=GRADIENT;
	    *---------------------------------------------------------*;
	    * End of Log-likelihood and Gradient function definitions *;
	    *---------------------------------------------------------*;
	    /*-------------------------------------------------------------*/
	    /*  Create Vector Containing Names of Parameters               */
	    /*-------------------------------------------------------------*/
	      Parameter = {"Intercept", %do J=1 %to &L; "&&VAR&J", %end;
	                   "Scale", "Shape (Kappa)" };
	         CNames = Parameter`;
	    /*-------------------------------------------------------------*/
	    /*                 Import the data into IML                    */
	    /*-------------------------------------------------------------*/
	     Use SUBJECT;                         /* one row per patient   */
	      read all var{&subject._r} into ID;     /* Subject ID numbers    */
	      read all var{&FIXED} into Z;        /* time-indep covariates */
	      read all var{K} into K;             /* obs per patient       */
	      read all var{&T1. &T2.} into T;     /* event/censor times    */
	      read all var{_CENSOR_} into CENSOR; /* censoring type        */
	      read all var{_START_ } into START ; /* start location of TDC */
	      read all var{_STOP_  } into STOP  ; /* stop  location of TDC */
	     Close SUBJECT;
	     *---------------------------------------------------*;
	     * Matrix TDC: multiple records per subject.         *;
	     *  Column 1 = Subject ID                            *;
	     *  Column 2 = Time when TDC changed value           *;
	     *  Column 3 to (3+Q) = TDC covariates               *;
	     *---------------------------------------------------*;
	     Use TDCVARS;
	      read all var{&subject._r &TDCTIME &TDC} into TDC;
	      read all var{&subject._r} into ID_DUPS;
	     Close TDCVARS;
	     *-------------------------------------------------------*;
	     * Set Initial Value of Parameter Vector (_Omega_)       *;
	     *-------------------------------------------------------*;
	     Use INITPARM;
	      read all var{ _SCALE_   } into _alpha_; _alpha_ = 1/_alpha_;
	      read all var{ Intercept } into _theta_; _theta_ = -_theta_#_alpha_;
	      read all var{ &FIXED    } into _gamma_; _gamma_ = -_gamma_#_alpha_;
	      read all var{ &TDC      } into _eta_  ; _eta_   = -_eta_#_alpha_  ;
	     Close INITPARM;
	       _Omega_ = ( _alpha_ || _theta_ || {0} || _gamma_ || _eta_ )`;
	     *-------------------------------------------------------------*;
	     *             Define a few useful quantities                  *;
	     *-------------------------------------------------------------*;
	      N = nrow(T);             /* number of patients               */
	      P = ncol(Z);             /* number of time-indep covariates  */
	      Q = ncol(TDC)-2;         /* number of time-dep covariates    */
	       r = P+Q+3;              /* Number of parameters             */
	      %IF %LENGTH(&FIXED_KAPPA.)^=0 %THEN %DO;
	      df = P + Q;              /* Full Model df, excluding scale,  */
	                               /* shape, and intercept parameters. */
	                               /* KAPPA is fixed.                  */
	      %END;
	      %ELSE %DO;
	      df = P + Q + 1 ;         /* Full Model df, excluding scale   */
	                               /* and intercept parameters.        */
	                               /* KAPPA is a free parameter to be  */
	                               /* estimated.                       */
	      %END;
	      NOBS = sum(K);           /* Total number of observations     */
	      EPSILON = abs(&EPSILON); /* Bound for Abs(1-KAPPA)           */

	    *---------------------------------------------------------*;
	    *  ML Estimation of Accelerated Failure Time TDC Model    *;
	    *---------------------------------------------------------*;
	    *----------------------------------------------------*;
	    *-  Set General Optimization options                -*;
	    *----------------------------------------------------*;
	      Options=J(1,11,.);
	      Options(|1|)=1;    /* 1=maximize, 0=minimize           */
	      %IF %UPCASE(&DETAILS)=SHORT %THEN %DO;
	        Options(|2|)=1;  /* Start, stop, iteration */
	      %END; %ELSE
	      %IF %UPCASE(&DETAILS)=LONG %THEN %DO;
	        Options(|2|)=2;  /* Start, stop, iteration, init & final parms */
	      %END; %ELSE
	      %DO;
	        Options(|2|)=0;  /* No optimization routine printout */
	      %END;
	      *--------------------------------------------------------------*;
	      *- Optimization options Specific to NLPFDD:finite differences -*;
	      *--------------------------------------------------------------*;
	      *--------------------------------------------------------------*;
	      *-  Option(|8|)= defines finite difference approximation used -*;
	      *-    to compute 1st or 2nd order derivatives                 -*;
	      *-             = ij                                           -*;
	      *-    where i=1 intervals by Gill, Murray, Saunders, & Wright -*;
	      *-              based on behavior of objective function.      -*;
	      *-          i=2 intervals based on behaviour of nonlinear     -*;
	      *-              constraint functions.                         -*;
	      *-          i=3 intervals based on behaviour of nonlinear     -*;
	      *-              constraint functions & objective function.    -*;
	      *-          j=0 fast but imprecise forward difference         -*;
	      *-         j^=0 expensive central difference formulas.        -*;
	      *--------------------------------------------------------------*;
	        Options(|8|)=11;
	      /*------------------------------------------------------------*/
	      /*     MLE of Parameter Vector is Omega0                      */
	      /*------------------------------------------------------------*/
	     %IF %LENGTH(&FIXED_KAPPA.)^=0 %THEN %DO;
	      /*------------------------------------------------------------*/
	      /*   MLE and Variances -- KAPPA is fixed and user-specified   */
	      /*------------------------------------------------------------*/
	      *--------------------------------------------------------------*;
	      *- Define constraint matrix.                                  -*;
	      *- Number of Columns = # parameters + 2 = r + 2 = (p+q+3) + 2 -*;
	      *- Number of Rows = #constraints + 2 = (p+q+2) + 2            -*;
	      *- Row #1 = Lower bound of r parameters (last 2 cols not used)-*;
	      *- Row #2 = Upper bound of r parameters (last 2 cols not used)-*;
	      *- Row p+q+2 = Constraints specified as follows:              -*;
	      *-          Columns 1:r = coefficient of parameters           -*;
	      *-          Column r+1 =  0 if relation is equality           -*;
	      *-                     =  1 if relation is >=                 -*;
	      *-                     = -1 if relation is <=                 -*;
	      *-          Column r+2 =  value of RHS of (in)equality        -*;
	      *--------------------------------------------------------------*;
	       *-------------------------------------------------------------------*;
	       *- Null Model (Intercept, Scale, and KAPPA Only) Constraint Matrix -*;
	       *-------------------------------------------------------------------*;
	  	  	%IF %LENGTH(&FIXED_SCALE.)^=0 %THEN %DO;
	       	*-------------------------------------------------------------------*;
	       	*- Set Options(|11|)= (p+q+2) equality constraints.                -*;
	       	*-   a) alpha=fixed                                                -*;
	       	*-   b) Kappa=fixed                                                -*;
	       	*-   c) P+Q parameters of covariates = 0                           -*;
	       	*-------------------------------------------------------------------*;
		       Options(|11|)=p+q+2;
			   /* alpha fixed, Kappa fixed, p+q parameters = 0 */;
			      Rows1_2 = J(2,1,&FIXED_SCALE.)|| J(2,1,.)||J(2,1,&FIXED_KAPPA.)||J(2,r-3,0)||J(2,2,.) ;
			   /* alpha = &FIXED_SCALE. */;
		          Row3    = {1} || J(1,r-1,0) || { 0 &FIXED_SCALE. } ;
			   /* p+q covariate parameters = 0, kappa = &FIXED_KAPPA. */;
			      Row_pq1 = J(p+q+1,2,0)||I(p+q+1)||J(p+q+1,1,0)||( {&FIXED_KAPPA.}//J(p+q,1,0) ) ;
			      Con = Rows1_2 // Row3 // Row_pq1;  
		  	%END;
		  	%ELSE %DO;
	       	*-------------------------------------------------------------------*;
	       	*- Set Options(|11|)= (p+q+1) equality + 1 inequality constraints. -*;
	       	*-   a) alpha > 0                                                  -*;
	       	*-   b) Kappa=fixed                                                -*;
	       	*-   c) P+Q parameters of covariates = 0                           -*;
	       	*-------------------------------------------------------------------*;
		       Options(|11|)=p+q+2;
			   /* alpha lower bound = &SCALE_LB, Kappa fixed, p+q parameters = 0 */;
			      Rows1_2 = {&SCALE_LB  . , .  .}||J(2,1,&FIXED_KAPPA.)||J(2,r-3,0)||J(2,2,.) ;
			   /* alpha >= 0 */;
		          Row3    = {1} || J(1,r-1,0) || { 1 0 } ;
			   /* p+q covariate parameters = 0, kappa = fixed */;
			      Row_pq1 = J(p+q+1,2,0)||I(p+q+1)||J(p+q+1,1,0)||( {&FIXED_KAPPA.}//J(p+q,1,0) ) ;
			      Con = Rows1_2 // Row3 // Row_pq1;  
		   	%END;
	       *------------------------------------------------------*;
	       *-  Initial values from Lifereg                       -*;
	       *------------------------------------------------------*;
	           init = _Omega_(|1:2|)//{&FIXED_KAPPA.}//J(r-3,1,0);
				%IF %LENGTH(&FIXED_SCALE.)^=0 %THEN %DO;
					init(|1|) = &FIXED_SCALE. ;
				%END;
	       *---------------------------------------------------------------------*;
	       *-  Estimate Null Model (Intercept, Scale, KAPPA (fixed)) Parameters -*;
	       *---------------------------------------------------------------------*;
	          call nlpnra(rc,omega0,"loglik",init,Options,con) grd="gradient";
	          call nlpfdd (F0,G0,H0,"loglik",omega0,Options,con) grd="gradient";
	 
	       *-----------------------------------------------------------------*;
	       *- Full Model(Kappa fixed, all other parameters to be estimated) -*;
	       *-----------------------------------------------------------------*;
	       *-  Initial values from Lifereg, KAPPA = user-specified value    -*;
	       *-----------------------------------------------------------------*;
	           init = _Omega_; init(|3|)=&FIXED_KAPPA.;
				%IF %LENGTH(&FIXED_SCALE.)^=0 %THEN %DO;
					init(|1|) = &FIXED_SCALE. ;
				%END;
	  	  	%IF %LENGTH(&FIXED_SCALE.)^=0 %THEN %DO;
				*------------------------------------------------*;
				*- Set Options(|11|)=2 to specify 2 constraints -*;
				*-     1) alpha = fixed value                   -*;
				*-     2) kappa = fixed value                   -*;
				*------------------------------------------------*;
				Options(|11|)=2;              
				Con=J(4,r+2,.);               /* 4 rows = 2 constraints + 2              */
				Con(|1,1|)  = &FIXED_SCALE. ; /* Lower  bound for alpha                  */
				Con(|2,1|)  = &FIXED_SCALE. ; /* Uppder bound for alpha                  */
				Con(|1,3|)  = &FIXED_KAPPA. ; /* Lower bound of 3rd parameter, Kappa     */
				Con(|2,3|)  = &FIXED_KAPPA. ; /* Upper bound of 3rd parameter, Kappa     */    
				Con(|3,1|)  = 1;    		  /* Coefficient of 1st parameter, alpha     */
				Con(|3,r+1|)= 0;  		      /* Col r+1 = 0 => relation is =            */
				Con(|3,r+2|)= &FIXED_SCALE.;  /* RHS of equality = &FIXED_SCALE          */ 
				Con(|4,3|)  = 1;    		  /* Coefficient of 3rd parameter, Kappa     */
				Con(|4,r+1|)= 0;  		      /* Col r+1 = 0 => relation is equality     */
				Con(|4,r+2|)= &FIXED_KAPPA.;  /* RHS of equality = &FIXED_KAPPA          */ 
	  	  	%END;
			%ELSE %DO;
				*------------------------------------------------*;
				*- Set Options(|11|)=2 to specify 2 constraints -*;
				*-     1) alpha > 0                             -*;
				*-     2) kappa = fixed value                   -*;
				*------------------------------------------------*;
				Options(|11|)=2;              
				Con=J(4,r+2,.);               /* 4 rows = 2 constraints + 2               */
				Con(|1,1|)= &SCALE_LB.    ;   /* Lower bound for alpha                    */
				Con(|1,3|)= &FIXED_KAPPA. ;   /* Lower bound of 3rd parameter, Kappa      */
				Con(|2,3|)= &FIXED_KAPPA. ;   /* Upper bound of 3rd parameter, Kappa      */    
				Con(|3,1|)= 1;    		      /* Coefficient of 1st parameter, alpha      */
				Con(|3,r+1|)= 1;  		      /* Col r+1 = 1 => relation is >=            */
				Con(|3,r+2|)= 0;              /* RHS of inequality = 0                    */ 
				Con(|4,3|)= 1;    		      /* Coefficient of 3rd parameter, Kappa      */
				Con(|4,r+1|)= 0;  		      /* Col r+1 = 0 => relation is equality      */
				Con(|4,r+2|)= &FIXED_KAPPA. ; /* RHS of equality = &FIXED_KAPPA           */ 
	  	  	%END;
	       *----------------------------------------------------------------------------------------------*;
	       *-  Estimate Full Model (Intercept, Scale, KAPPA (fixed), Fixed and TD covariates) Parameters -*;
	       *----------------------------------------------------------------------------------------------*;
	          call nlpnra(rc,omega,"loglik",init,Options,con) grd="gradient";
	          call nlpfdd (F,G,H,"loglik",omega,Options) grd="gradient";
	          *--------------------*;
	          * Information Matrix *;
	          *--------------------*;
	             Info=-H;
	     %END;
	     %ELSE %DO;
	      /*-------------------------------------------------------------*/
	      /*   MLE and Variances -- KAPPA to be estimated from the model */
	      /*-------------------------------------------------------------*/
	      *--------------------------------------------------------------*;
	      *- Define constraint matrix.                                  -*;
	      *- Number of Columns = # parameters + 2 = r + 2 = (p+q+3) + 2 -*;
	      *- Number of Rows = #constraints + 2 = (p+q+2) + 2            -*;
	      *- Row #1 = Lower bound of r parameters (last 2 cols not used)-*;
	      *- Row #2 = Upper bound of r parameters (last 2 cols not used)-*;
	      *- Row p+q+2 = Constraints specified as follows:              -*;
	      *-          Columns 1:r = coefficient of parameters           -*;
	      *-          Column r+1 =  0 if relation is equality           -*;
	      *-                     =  1 if relation is >=                 -*;
	      *-                     = -1 if relation is <=                 -*;
	      *-          Column r+2 =  value of RHS of (in)equality        -*;
	      *--------------------------------------------------------------*;
	       *----------------------------------------------------------------------*; 
	       *- Null Model (Intercept, Scale, and KAPPA_Ho Only) Constraint Matrix -*;
	       *----------------------------------------------------------------------*; 
	       *-------------------------------------------------------------------*;
	       *- Set Options(|11|)= (p+q+1) equality + 1 inequality constraints. -*;
	       *-   a) alpha > 0                                                  -*;
	       *-   b) Kappa = Kappa_Ho                                           -*;
	       *-   c) P+Q parameters of covariates = 0                           -*;
	       *-------------------------------------------------------------------*;
	       Options(|11|)=p+q+2;
		   /* alpha lower bound = 0, Kappa = Kappa_Ho, p+q parameters = 0 */;
		      Rows1_2 = {&SCALE_LB  . , .  .}||J(2,1,&KAPPA_Ho.)||J(2,r-3,0)||J(2,2,.) ;
		   /* alpha >= 0 */;
	          Row3    = {1} || J(1,r-1,0) || { 1 0 } ;
		   /* p+q covariate parameters = 0, kappa = Kappa_Ho */;
		      Row_pq1 = J(p+q+1,2,0)||I(p+q+1)||J(p+q+1,1,0)||( {&KAPPA_Ho.}//J(p+q,1,0) ) ;
	       Con = Rows1_2 // Row3 // Row_pq1;  
	       *----------------------------------*;
	       *-  Initial values from Lifereg   -*;
	       *----------------------------------*;
	           init = _Omega_(|1:2|)//{&KAPPA_Ho.}//J(r-3,1,0);
	       *------------------------------------------------------------------------*;
	       *-  Estimate Null Model (Intercept, Scale, KAPPA = Kappa_Ho) Parameters -*;
	       *------------------------------------------------------------------------*;
	          call nlpnra(rc,omega0,"loglik",init,Options,con) grd="gradient";
	          call nlpfdd (F0,G0,H0,"loglik",omega0,Options,con) grd="gradient";
	       *-----------------------------------------------------------------*;
	       *- Full Model(All parameters to be estimated)                    -*;
	       *-----------------------------------------------------------------*;
	       *--------------------------------------------------------------------------*;
	       *-  Initial values from Lifereg, initial value of KAPPA is user-specified -*;
	       *--------------------------------------------------------------------------*;
	           init = _Omega_; init(|3|)=&INIT_KAPPA.;
	       *------------------------------------------------*;
	       *- Set Options(|11|)=3 to specify 3 constraints -*;
		   *-     1) alpha > 0                             -*;
		   *-     2) kappa >= KAPPA_LB                     -*;
		   *-     2) kappa <= KAPPA_UB                     -*;
	       *------------------------------------------------*;
			Options(|11|)=3;              
			Con=J(5,r+2,.);               /* 5 rows = 3 constraints + 2               */
			Con(|1,1|)= &SCALE_LB.    ;   /* Lower bound for alpha                    */
			Con(|1,3|)= &KAPPA_LB.    ;   /* Lower bound of 3rd parameter, Kappa      */
			Con(|2,3|)= &KAPPA_UB.    ;   /* Upper bound of 3rd parameter, Kappa      */    
			Con(|3,1|)  =  1;   	      /* Coefficient of 1st parameter, alpha      */
			Con(|3,r+1|)=  1;  		      /* Col r+1 = 1 => relation is >=            */
			Con(|3,r+2|)=  0;             /* RHS of inequality = 0                    */ 
			Con(|4,3|)  =  1;    	      /* Coefficient of 3rd parameter, Kappa      */
			Con(|4,r+1|)=  1;  		      /* Col r+1 = 1 => relation is >=            */
			Con(|4,r+2|)=  &KAPPA_LB. ;   /* RHS of equality = &KAPPA_LB              */ 
			Con(|5,3|)  =  1;  		      /* Coefficient of 3rd parameter, Kappa      */
			Con(|5,r+1|)= -1;  		      /* Col r+1 = -1 => relation is <=           */
			Con(|5,r+2|)= &KAPPA_UB.  ;   /* RHS of equality = &KAPPA_UB              */ 
	       *--------------------------------------------------------------------------------------*;
	       *-  Estimate Full Model (Intercept, Scale, KAPPA, Fixed and TD covariates) Parameters -*;
	       *--------------------------------------------------------------------------------------*;
	           call nlpnra(rc,omega,"loglik",init,Options,Con) grd="gradient";
	           call nlpfdd (F,G,H,"loglik",omega,Options,Con) grd="gradient";
	          *--------------------*;
	          * Information Matrix *;
	          *--------------------*;
	             Info = -H;
	     %END;
	      /*------------------------------------------------------------*/
	      /* Rearrange Order of Parameter Vector and Covariance Matrix  */
	      /*------------------------------------------------------------*/
	          *----------------------*;
	          * Move Alpha to bottom *;
	          *----------------------*;
	           PostMult=(J(1,r-1,0)||{1}   ) // (I(r-1)||J(r-1,1,0));
	            PreMult=(J(r-1,1,0)||I(r-1)) // (   {1}||J(1,r-1,0));
	              Parms=Omega*Postmult;
				  Info=Premult*Info*Postmult;
	          *-----------------------------------*;
	          * Move Kappa to bottom, after Alpha *;
	          *-----------------------------------*;
	           PostMult=({1}//J(r-1,1,0))||(J(2,r-2,0)//I(r-2))||
	                    ({0,1}//J(r-2,1,0));
	            PreMult=({1}//J(r-1,1,0))||(J(r-1,1,0)//{1})||
	                    (J(1,r-2,0)//I(r-2)//J(1,r-2,0));
	              Parms=(Parms*Postmult)`;
				  Info=Premult*Info*Postmult;
	      /*------------------------------------------------------------*/
	      /* Compute Model LRT on P+Q degress-of-freedom.               */
	      /*------------------------------------------------------------*/
	             LogLNull = F0; LogL=F;
	               LRT_X2 = -2#(F0-F);
	                LRT_P = 1-Probchi(LRT_X2,df);
	      /*------------------------------------------------------------*/
	      /* Compute Wald p-values and 95% CI of Parameter Estimates    */
	      /*------------------------------------------------------------*/
			%IF %LENGTH(&FIXED_KAPPA.)^=0 %THEN %DO;
				%IF %LENGTH(&FIXED_SCALE.)^=0 %THEN %DO;
	             SE_Parms = SE_parms // {0,0};
	              Low95CI = Low95CI  // {.,.};
			     High95CI = High95CI // {.,.};
	              Wald_X2 = Wald_X2  // {.,.};
	               Wald_P = Wald_P   // {.,.};
				%END;
				%ELSE %DO;
	             SE_Parms = SE_parms // {0};
	              Low95CI = Low95CI  // {.};
			     High95CI = High95CI // {.};
	              Wald_X2 = Wald_X2  // {.};
	               Wald_P = Wald_P   // {.};
				%END;
			%END;
			%ELSE %DO;
			  Parms_Ho = J(r-1,1,0)//{&kappa_Ho.};
			%END;
	      /*-----------------------*/
	      /* Create Output Datsets */
	      /*-----------------------*/
	      Estimate=Parms; 
		  %if &nreps=0 %then %do;
		      Create &SAVELIB..LRTMODEL_R&nreps VAR {N NObs LogLNull LogL df
		                                     LRT_X2 LRT_P };
		         append VAR {N NObs LogLNull LogL df LRT_X2 LRT_P };
		         close &SAVELIB..LRTMODEL_R&nreps;
		  %end;

	      Create &SAVELIB..PARMEST_R&nreps
	         VAR {Parameter Estimate };
	      append VAR {Parameter Estimate };
	         close &SAVELIB..PARMEST_R&nreps;
	   QUIT;
		proc datasets library=work memtype=data nolist;
			delete InputDat ReducedData Check_T0 censored_int censored_oth dummy1 dummy2 Recs_Added
				   UseData newuse TDCVARS SUBJECT INITPARM ;
		quit;
	%END; /*End estimation of resampled data*/

   *----------------------------------------------------------*;
   * End IML Estimation of Accelerated Failure Time-TDC Model *;
   *----------------------------------------------------------*;

   *------------------------------*;
   * Collate Results for Printing *;
   *------------------------------*;
    proc format;
       PICTURE  PFMT      .=' '
               LOW-<0.0001 ='<0.0001' (NOEDIT)
            0.0001-HIGH    =' 9.9999';
       PICTURE  X2FMT     .=' '
                 0-HIGH    ='0,009.9999';
       VALUE    LINEFMT   1="__________"
                          2="___________"
                          3="______________________"
                          4="__________________";
       PICTURE  CILOFMT .=' ' LOW-<0='0,0009.9999 ,' (PREFIX='( -')
                              0-HIGH='0,009.9999 ,'  (PREFIX='(  ');
       PICTURE  CIHIFMT .=' ' LOW-<0='0,009.9999 )'   (PREFIX='-' )
                              0-HIGH='0,009.9999 )'   (PREFIX=' ' );
    run;
    data _LRT;
       set &SAVELIB..LRTMODEL_R0; 
    run;
    data _PARMEST_BS;
       set PARMEST_R1-PARMEST_R&reps;
    run;
	proc means data=_PARMEST_BS noprint;
		class parameter;
		types parameter;
		var estimate;
		output out=stat_bydata std=SE ;
	run;

	/*construct new 95% with SE from bootstrapping*/
	%sortby(&SAVELIB..PARMEST_R0, parameter);
	%sortby(stat_bydata, parameter);
	data _PARMEST;
		merge &SAVELIB..PARMEST_R0 stat_bydata;
		by parameter;

		Low95CI = ROUND(estimate-1.96*SE,0.0001);
		High95CI  = ROUND(estimate+1.96*SE,0.0001);

		array stat estimate Low95CI High95CI;
		array expv relrisk exp_Low95CI exp_High95CI;
		do over expv;
			expv = exp(stat);
		end;
	run;



	data ds_out;
		set _LRT _PARMEST;
	run;
   *-----------------*;
   * Display Results *;
   *-----------------*;
	options orientation=portrait topmargin=0.5in bottommargin=0.5in leftmargin=0.5in rightmargin=0.5in
			ls=112 ps=60;
	data _null_;
	  set ds_out;
	  file print header=newpage;
	  if _n_=1 and LogLNull =. then do;
	        put //		@43 "Kappa" @57 "="
							%if %length(&fixed_kappa.)^= 0 %then %do;
									@60 "&fixed_kappa. (user specified)"
								%if %length(&fixed_scale.)^= 0 %then %do;
								/	@43 "Scale" @57 "=" @60 "&fixed_scale. (user specified)"
								%end;
								/	@43 "Distribution" @57 "=" @60 "(No name)"
									%if &fixed_kappa. = 1 %then @60 "Log-logistic";
									%if &fixed_kappa. = 0 %then %do;
											%if %length(&fixed_scale.) = 0 %then @60 "Weibull     ";
											%if %length(&fixed_scale.)^= 0 %then @60 "Exponential ";
									%end;
							%end;
							%if %length(&fixed_kappa.) = 0 %then %do;
									@60 "(unspecified)"
								/	@43 "Distribution" @57 "=" 
							%end;
				//		@40 "****************************************"
				/ 		@40 "** MODEL ESTIMATION DID NOT CONVERGE. **"
				/		@40 "**  TRY A DIFFERENT VALUE OF KAPPA.   **"
				/		@40 "****************************************"
			;
			abort;
	  end;
	  else do;
		if _n_=1 then do;
		%if %length(&fixed_kappa.)^= 0 %then %do;
	        put //		@43 "Kappa" @57 "="
									@60 "&fixed_kappa. (user specified)"
								%if %length(&fixed_scale.)^= 0 %then %do;
								/	@43 "Scale" @57 "=" @60 "&fixed_scale. (user specified)"
								%end;
								/	@43 "Distribution" @57 "=" @60 "(No name)"
									%if &fixed_kappa. = 1 %then @60 "Log-logistic";
									%if &fixed_kappa. = 0 %then %do;
											%if %length(&fixed_scale.) = 0 %then @60 "Weibull     ";
											%if %length(&fixed_scale.)^= 0 %then @60 "Exponential ";
									%end;
		%end;
		%else %do;
	        put //		@43 "Kappa (Null)" @57 "="
									@60 "&kappa_Ho." 
								/	@43 "Distribution" @57 "=" @60 "(No name)"
									%if &kappa_ho. = 1 %then   @60 "Log-logistic";
									%if &kappa_ho. = 0 %then   @60 "Weibull     ";
		%end;
 
				///		@49 "Log-Likelihood"
						@74 "Likelihood Ratio Test"
			overprint	@43 27*"_" @71 29*"_"
					//	@15 "Subjects (N)"
						@29 "Observations"
						@45 "Null Model"
						@59 "Full Model"
						@73 "df"
/*						@79 "Chi-Square"*/
/*						@91 "p-value"*/
					//	@15 N			comma12.0
						@29 NObs		comma12.0
						@43 LogLNull	best.
						@57 LogL		best.
						@71 DF			4.0
						@77 LRT_X2		best.
						@91LRT_P		pfmt.
 
 
				////	@43 "Maximum Likelihood Estimates"
			overprint	@43 28*"_" 
 
				//		@55 "95% Confidence Interval"
			overprint	@53 27*"_" 
				//		@9  "Parameter"
						@29 "Estimate"
						@41 "Std Err"
						@54 "Lower Limit"
						@68 "Upper Limit"
				//
		;
		end;
		else do;
			put			@9  Parameter
						@27 Estimate	comma10.4
						@40 SE			comma8.4
						@50 Low95CI		cilofmt.
						@67 High95CI	cihifmt.
			;
		end;
		return;
	  end;
	  newpage:
        put  // @13 "*" +1 85*"-" @100 "*"
              / @13 "*" @15 "Case-cohort interval-censored data with time-dependent covariates"
						@100 "*"
			  / @13 "*" @20 "by Gao X., Hudgens M., and Zou, F. - (2021+)"
						@100 "*"
			  / @13 "*" @25 "CCOH_PARM_ICE (Version 1.0): SAS Macro developed by Xiaoming Gao "
			  / @13 "*" @25 "based on PARM_ICE.SAS (Version 3.2) written by Oliver M. Bautista"
						@100 "*"
			  / @13 "*" +1 85*"-" @100 "*"
			;
		return;
	run;
   *-------------------------------------------------------------*;
   * Erase temporary SAS datasets created in the SASWORK library *;
   *-------------------------------------------------------------*;
	proc datasets library=work memtype=data nolist;
		delete InputDat ReducedData Check_T0 censored_int censored_oth dummy1 dummy2 Recs_Added
			   UseData newuse TDCVARS SUBJECT INITPARM _LRT  ds_out
				rsample widersample interval_bs interval_rp interval_hits interval_hits2 ds_id
				LRTMODEL_R0 PARMEST_R0-PARMEST_R&reps _PARMEST_BS stat_bydata;
	quit;
%MEND;
 
