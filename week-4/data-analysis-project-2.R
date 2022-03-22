################################################################################
# Data analysis project - carry out an IPTW causal analysis
################################################################################

# For this assignment we will use data from Lalonde (1986),
# that aimed to evaluate the impact of  National
# Supported Work (NSW) Demonstration, which is a labor training program, 
# on post-intervention income levels. Interest is in estimating the causal 
# effect of this training program on income.

# packages
library(tableone)
library(Matching)
library(ipw)
library(survey)
library(tidyverse)

# load data
library(MatchIt)
data(lalonde)
lalonde <- tibble(lalonde)

lalonde <- lalonde %>% 
  mutate(black = if_else(race == "black", 1, 0),
         hispan = if_else(race == "hispan", 1, 0)) %>% 
  dplyr::select(-race)


# The data have n=614 subjects and 10 variables
# 
# `age` age in years. 
# `educ` years of schooling. 
# `black` indicator variable for blacks. 
# `hispan` indicator variable for Hispanics. 
# `married` indicator variable for marital status. 
# `nodegree` indicator variable for high school diploma. 
# `re74` real earnings in 1974. 
# `re75` real earnings in 1975. 
# `re78` real earnings in 1978. 
# `treat` an indicator variable for treatment status.
#
# The outcome is
# re78 – post-intervention income.
# 
# The treatment is
# treat – which is equal to 1 if the subject received the 
# labor training and equal to 0 otherwise.
# 
# The potential confounding
# variables are: age, educ, black, hispan, married, nodegree, re74, re75.





# Question 1 -------------------------------------------------------------------

# Fit a propensity score model. Use a logistic regression model, where the 
# outcome is treatment. Include the 8 confounding variables in the model as 
# predictors, with no interaction terms or non-linear terms (such as squared 
# terms). Obtain the propensity score for each subject. Next, obtain the
# inverse probability of treatment weights for each subject.

# What are the minimum and maximum weights? 

# fit a logistic model
model <- glm(treat ~ 
               age +
               educ + 
               black + 
               hispan + 
               married + 
               nodegree + 
               re74 +
               re75,
             data = lalonde,
             family = binomial(link = "logit"))

#show coefficients etc
summary(model)

## value of propensity score for each subject
ps <- predict(model, type = "response")

#create weights
weight <- ifelse(lalonde$treat == 1, 1/(ps), 1/(1-ps))


min(weight)
max(weight)


###################
# 1.01 and 40.1   #
###################





# Question 2 -------------------------------------------------------------------
# Find the standardized differences for each confounder on the weighted 
# (pseudo) population. What is the standardized difference for nodegree? 

# confounding variables
xvars <- c("age", "educ", "black", "hispan", "married", "nodegree","re74", "re75")

#apply weights to data
weighteddata <- survey::svydesign(ids = ~ 1, data = lalonde, weights = ~ weight)


#weighted table
weightedtable <- svyCreateTableOne(vars = xvars, 
                                   strata = "treat", 
                                   data = weighteddata, 
                                   test = FALSE)
## Show table with SMD
print(weightedtable, smd = TRUE)


###################
# 0.11            #
###################










# Question 3 -------------------------------------------------------------------
# Using IPTW, find the estimate and 95% confidence interval for the average 
# causal effect. This can be obtained from svyglm

#fit a marginal structural model (risk difference)
msm <- (svyglm(re78 ~ treat, 
               design = svydesign(~ 1, weights = ~weight,
                                  data = lalonde)))
coef(msm)
confint(msm)


############################################
# Est: 224.68 95% CI: (-1559.32, 2008.67)  #
############################################









# Question 4 -------------------------------------------------------------------
# Now truncate the weights at the 1st and 99th percentiles. This can be done
# with the trunc=0.01 option in svyglm.

# Using IPTW with the truncated weights, find the estimate and 95% confidence 
# interval for the average causal effect


# fit propensity score model to get weights, but truncated
weightmodel <- ipwpoint(exposure = treat, 
                        family = "binomial", 
                        link ="logit",
                        denominator =~ 
                          age +      
                          educ +     
                          black +   
                          hispan +   
                          married +  
                          nodegree + 
                          re74 +     
                          re75, 
                        trunc=.01,
                        data = lalonde %>% as.data.frame) #ipw cannot handle tibbles



lalonde$wt <- weightmodel$weights.trun

#fit a marginal structural model (risk difference)
msm <- (svyglm(re78 ~ treat, 
               design = svydesign(~ 1, weights = ~wt,
                                  trunc=0.01,
                                  data = lalonde)))
coef(msm)
confint(msm)



############################################
# 486.93  (-1090.64, 2064.51)              #
############################################
 


















