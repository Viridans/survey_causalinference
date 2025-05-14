/* IPTW with survey weights in SAS

We are interested in examining the effect of healthcare access documented in the 2013-2014
waves of NHANES on individual mortality, assessed in 2019 from NDI data. Confounders
include gender, age, race/ethnicity, education and income.

Variables

Exposure: healthcare_access
Outcome: mortstat
Covariates: gender, age_cat, race_ethnicity, education, income

We will use IPW to get an average treatment effect (ATE) estimate. To generalize this
estimate to the population that CDC sampled from, we'll incorporate information about
the hierarchical sampling process:

Strata: strata
Cluster: psu
Individual: seqn
Individual weight: svyw

*/

/* Import data */

PROC IMPORT out=nhanes
datafile="yourpathhere\nhanes.csv"
	DBMS=CSV REPLACE;
    GETNAMES=YES;
    DATAROW=2; 
RUN;

/* Examine data (we're using a complete case dataset here for example purposes, N=5480)*/

PROC FREQ data=nhanes; *variable values and proportions;
	TABLES healthcare_access mortstat gender age_cat race_ethnicity education income /missing;
RUN;

PROC SQL; *it's always a good idea to check for duplicate individuals;
	create table counts as
	select seqn, count(*) as IDCount from nhanes
	group by seqn;

PROC PRINT data=counts;
	WHERE IDCount > 1;
RUN;

/* Assuming we want to make inference about our source population instead of our sample
(generally what we're trying to do when answering public health questions), there are
two schools of thought about how to handle the treatment (PS) model:


- One states that propensity scores pertain to the sample's covariate distribution,
not the source population's--they allow for calculation of a marginal estimate
among the OBSERVED treated and untreated. Therefore, survey weights should be a
covariate in the propensity score calculation, and an unweighted model (e.g.
logistic regression) should be used to calculate propensity scores.

- Another states that since weight generation is a process that occurs after
treatment (exposure) occurs, and there's the possibility that treatment status
is associated with selection into the study, we need to weight back to the source
population first through using a weighted model to calculate propensity scores--
essentially work backward.

To examine this, we'll do four models, where the last two reflect these options:


Model 1: IPW without any inclusion of the survey structure in the treatment or outcome
models, for comparison.

Model 2: No inclusion of weights in the treatment model; weighted outcome model.

Model 3: Generation of IPW weights with the survey weights as a covariate in the treatment
model; weighted outcome model.

Model 4: Generation of IPW weights using a weighted regression in the treatment model;
weighted outcome model.

*/


/***********************************

			 Model 1

***********************************/


/* First, generate the propensity scores. We can do this using PROC PSMATCH, which
uses logistic regression (see: "Applying Propensity Score Methods to Complex Survey Data"
by Patrick Karabon: https://support.sas.com/resources/papers/proceedings19/3634-2019.pdf)

In the output dataset [getate], variable [_ps_] is the calculated propensity score
and variable [_ate_] is the ipw weight, calculated in relation to the entire
dataset (ATT would be in relation to the treated.) To get weights for the ATT
instead, you would specify the "ATTWGT" option.

*/

PROC PSMATCH data=nhanes region=allobs;
	CLASS healthcare_access gender age_cat race_ethnicity education income;
	PSMODEL healthcare_access(treated="1") = gender age_cat race_ethnicity education income;
	OUTPUT OUT(obs=all)=getate ATEWGT=_ATE_;
RUN;

/* Now, run PROC GENMOD to get an estimate for the sample, applying IPW weights.

(REPEATED SUBJECT allows for calculation of the correct SEs.) */

PROC GENMOD data=getate DESCENDING;
	CLASS seqn;
	MODEL mortstat = healthcare_access /link=log distribution=binomial;
	WEIGHT _ATE_;
	REPEATED SUBJECT = seqn;
	ESTIMATE 'ATE RR' healthcare_access 1;
RUN;


/***********************************

			 Model 2

***********************************/

/* Here, we'll incorporate the survey structure into the outcome model only. To do this,
we will combine the IPTW weights calculated above and the survey weights through multiplying,
the way we'd do it when doing IPTW for confounder control and IPTW for missingness. */

DATA final;
	MERGE nhanes getate;
	BY seqn;
	finalwt = _ATE_*svyw;
RUN;

/* SAS has built-in procedures that incorporate survey variables for ordinary regression
(PROC SURVEYREG) and for logistic regression (PROC SURVEYLOGISTIC) but as of v.9.4
does not have a survey version of GENMOD. Here, I used a macro developed by Alan
Ricardo da Silva and extended by Noah Padgett, which yields a similar result
to R routine svyglm.

Original conference paper: https://support.sas.com/resources/papers/proceedings17/0268-2017.pdf
Github link: https://github.com/noah-padgett/surveygenmod2

(There is one small error in the .sas file--make sure the %mend statement reads
%mend surveygenmod2;)

The outcome is binary in this case, so I've used a binomial distribution with log link.
Supplied variables must be numeric, and if you're not doing a complete case analysis,
missing values should be imputed. y= is outcome, x= is exposure, weight = the previously
calculated final weight, and other sampling variables go into strata= and cluster=.

*/

%include 'yourpathhere\surveygenmod2.sas';

%surveygenmod2(data=final, y=mortstat, x=healthcare_access, domain=, dist=binomial, link=log
,	weight=finalwt, strata=strata, cluster=psu, fpc=, scale=, intercept=y, vadjust=n
,	output=paramEST, output2=vcovEst);


/***********************************

			 Model 3

***********************************/

/* This approach adds the survey weights as a covariate in the treatment
model and then runs a weighted regression outcome model, as above. */

PROC PSMATCH data=nhanes region=allobs;
	CLASS healthcare_access gender age_cat race_ethnicity education income;
	PSMODEL healthcare_access(treated="1") = gender age_cat race_ethnicity education income svyw; *notice svyw;
	OUTPUT OUT(obs=all)=getate ATEWGT=_ATE_;
RUN;

DATA final; *multiply ipw and survey weights;
	MERGE nhanes getate;
	BY seqn;
	finalwt = _ATE_*svyw;
RUN;

%surveygenmod2(data=final, y=mortstat, x=healthcare_access, domain=, dist=binomial, link=log
,	weight=finalwt, strata=strata, cluster=psu, fpc=, scale=, intercept=y, vadjust=n
,	output=paramEST, output2=vcovEst); *outcome model;


/***********************************

			 Model 4

***********************************/

/* In this approach, we calculate the PS scores using a weighted regression (PROC SURVEYLOGISTIC),
instead of applying the weight as a covariate, then proceed with a weighted regression outcome model. */

PROC SURVEYLOGISTIC data=nhanes; *get PS scores;
	STRATA strata;
	CLUSTER psu;
	WEIGHT svyw;
	CLASS healthcare_access gender age_cat race_ethnicity education income;
	MODEL healthcare_access(event='1') = gender age_cat race_ethnicity education income;
	OUTPUT OUT=weightedps predicted=PS;
RUN;

DATA ipw; *get ipw weights;
	SET weightedps;
	IF healthcare_access = 1 THEN ipw = 1/ps;
	ELSE IF healthcare_access = 0 THEN ipw = 1/(1-ps);
RUN;

DATA final; *multiply ipw and survey weights;
	MERGE nhanes ipw;
	BY seqn;
	finalwt = ipw*svyw;
RUN;

%surveygenmod2(data=final, y=mortstat, x=healthcare_access, domain=, dist=binomial, link=log
,	weight=finalwt, strata=strata, cluster=psu, fpc=, scale=, intercept=y, vadjust=n
,	output=paramEST, output2=vcovEst); *outcome model;
