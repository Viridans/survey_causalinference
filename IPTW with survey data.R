# IPTW with survey weights in R
# 
# We are interested in examining the effect of healthcare access documented in the 2013-2014
# waves of NHANES on individual mortality, assessed in 2019 from NDI data. Confounders
# include gender, age, race/ethnicity, education and income.
# 
# Variables
# 
# Exposure: healthcare_access
# Outcome: mortstat
# Covariates: gender, age_cat, race_ethnicity, education, income
# 
# We will use IPW to get an average treatment effect (ATE) estimate. To generalize this
# estimate to the population that CDC sampled from, we'll incorporate information about
# the hierarchical sampling process:
# 
# Strata: strata
# Cluster: psu
# Individual: seqn
# Individual weight: svyw


# Import data and load packages

library(tidyverse)
library(survey)

nhanes <- read.csv("yourpathhere/nhanes.csv")


# Assuming we want to make inference about our source population instead of our sample
# (generally what we're trying to do when answering public health questions), there are
# two schools of thought about how to handle the treatment (PS) model:
# 
# 
# - One states that propensity scores pertain to the sample's covariate distribution,
#   not the source population's--they allow for calculation of a marginal estimate
#   among the OBSERVED treated and untreated. Therefore, survey weights should be a
#   covariate in the propensity score calculation, and an unweighted model (e.g.
#   logistic regression) should be used to calculate propensity scores.
# 
# - Another states that since weight generation is a process that occurs after
#   treatment (exposure) occurs, and there's the possibility that treatment status
#   is associated with selection into the study, we need to weight back to the source
#   population first through using a weighted model to calculate propensity scores--
#   essentially work backward.
#   
# To examine this, we'll do four models, where the last two reflect these options:
# 
# 
# Model 1: IPW without any inclusion of the survey structure in the treatment or outcome
# models, for comparison.
# 
# Model 2: No inclusion of weights in the treatment model; weighted outcome model.
# 
# Model 3: Generation of IPW weights with the survey weights as a covariate in the treatment
# model; weighted outcome model.
# 
# Model 4: Generation of IPW weights using a weighted regression in the treatment model;
# weighted outcome model.


##########################
#                        #
#        Model 1         #
#                        #
##########################

# First, generate the propensity scores.
# We use the glm() function with the family="binomial" option to request
# logistic regression.

# Treatment model
ps_fit <- glm(healthcare_access ~ gender +
                                  age_cat +
                                  race_ethnicity +
                                  education +
                                  income,
              data = nhanes,
              family = "binomial")

# Get propensity scores
nhanes$ps <- predict(ps_fit, type="response")

ipw <- nhanes %>%
  mutate(ipw = NA) %>%
  mutate(ipw = ifelse(healthcare_access == 1, 1/ps, 1/(1 - ps)))

# We use the survey functions to incorporate the IPTW weight only.
# "id" refers to the cluster variable, "strata" to the strata variable, and
# "weights" to the individual weights.
# Set "id" to 1, "strata" to NULL, and specify the IPTW weight for "weights".

# Specify design variables
effect_nhanes <- svydesign(id = ~1, 
                           strata  = NULL,
                           weights = ~ipw,
                           data    = ipw)

# Call outcome model.
# Use the "quasibinomial(link=log)" option to request log-binomial regression.

effect <- svyglm(mortstat ~ healthcare_access, 
                 design = effect_nhanes, 
                 family = quasibinomial(link = log))

exp(coefficients(effect))
exp(confint.default(effect))


##########################
#                        #
#        Model 2         #
#                        #
##########################

# Here, we'll incorporate the survey structure into the outcome model only. To do this,
# we will combine the IPTW weights calculated above and the survey weights through multiplying,
# the way we'd do it when doing IPTW for confounder control and IPTW for missingness.

ipw <- nhanes %>%
  mutate(ipw = NA) %>%
  mutate(ipw = ifelse(healthcare_access == 1, 1/ps, 1/(1 - ps))) %>%
  mutate(svyw_ipw = svyw*ipw)

# We use the survey functions to incorporate the sampling design variables.
# "id" refers to the cluster variable, "strata" to the strata variable, and
# "weights" to the individual weights.
# The "nest=TRUE" option tells R that the clusters are nested within strata.

# Specify design variables
survey_data <- svydesign(id = ~psu,
                         strata  = ~strata,
                         weights = ~svyw_ipw,
                         nest    = TRUE,
                         data    = ipw)

# Now run a weighted regression for the outcome model, including the design
# object "survey_data".
# Use the "quasibinomial(link=log)" option to request log-binomial regression.

# Outcome model
fit <- svyglm(mortstat ~ healthcare_access, 
                 design = survey_data, 
                 family = quasibinomial(link=log))
exp(coefficients(fit))
exp(confint.default(fit))


##########################
#                        #
#        Model 3         #
#                        #
##########################

# This approach adds the survey weights as a covariate in the treatment
# model and then runs a weighted regression outcome model, as above.

# Treatment model
ps_fit <- glm(healthcare_access ~ gender +
                age_cat +
                race_ethnicity +
                education +
                income +
                svyw,
              data = nhanes,
              family = "binomial")

# Get propensity scores
nhanes$ps <- predict(ps_fit, type="response")

# Get IPTW weights and multiply to svyw
ipw <- nhanes %>%
  mutate(ipw = NA) %>%
  mutate(ipw = ifelse(healthcare_access == 1, 1/ps, 1/(1 - ps))) %>%
  mutate(svyw_ipw = svyw*ipw)

# Outcome model
survey_data <- svydesign(id = ~psu,
                         strata  = ~strata,
                         weights = ~svyw_ipw,
                         nest    = TRUE,
                         data    = ipw)

fit <- svyglm(mortstat ~ healthcare_access, 
              design = survey_data, 
              family = quasibinomial(link=log))
exp(coefficients(fit))
exp(confint.default(fit))


##########################
#                        #
#        Model 4         #
#                        #
##########################

# In this approach, we calculate the PS scores using a weighted regression (svyglm)
# instead of applying the weight as a covariate, then proceed with a weighted regression
# outcome model.
# Do not specify "(link=log)", only "quasibinomial" alone, to use weighted logistic
# regression in the PS calculation.

# Treatment model
survey_data <- svydesign(id = ~psu,
                         strata  = ~strata,
                         weights = ~svyw,
                         nest    = TRUE,
                         data    = nhanes)

ps_fit <- svyglm(healthcare_access ~ gender +
                                  age_cat +
                                  race_ethnicity +
                                  education +
                                  income,
              design = survey_data, 
              family = quasibinomial)

# Get propensity scores
nhanes$ps <- predict(ps_fit, type="response")

# Get IPTW weights and multiply to svyw
ipw <- nhanes %>%
  mutate(ipw = NA) %>%
  mutate(ipw = ifelse(healthcare_access == 1, 1/ps, 1/(1 - ps))) %>%
  mutate(svyw_ipw = svyw*ipw)

# Outcome model
survey_data <- svydesign(id = ~psu,
                         strata  = ~strata,
                         weights = ~svyw_ipw,
                         nest    = TRUE,
                         data    = ipw)

fit <- svyglm(mortstat ~ healthcare_access, 
              design = survey_data, 
              family = quasibinomial(link=log))
exp(coefficients(fit))
exp(confint.default(fit))
