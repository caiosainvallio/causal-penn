################################################################################
# Data analysis project - analyze data in R using propensity score matching
################################################################################



# For this assignment we will use data from Lalonde (1986),
# that aimed to evaluate the impact of  National
# Supported Work (NSW) Demonstration, which is a labor training program, on
# post-intervention income levels. Interest is in estimating the causal effect of
# this training program on income.


# packages
library(tableone)
library(Matching)
library(MatchIt)
library(tidyverse)

# load data
data("lalonde")
lalonde <- tibble::tibble(lalonde)
lalonde

lalonde <- lalonde %>% 
  mutate(black = if_else(race == "black", 1, 0),
         hispan = if_else(race == "hispan", 1, 0)) %>% 
  select(-c(race))


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
# Find the standardized differences for all of the confounding variables 
# (pre-matching). What is the standardized difference for married (to nearest hundredth)?

# confounding variables
xvars <- c("age", "educ", "black", "hispan", "married", "nodegree","re74", "re75")

#look at a table
table <- tableone::CreateTableOne(vars = xvars,
                                   strata = "treat", 
                                   data = lalonde,
                                   test = T)

## include standardized mean difference (SMD)
print(table, smd=TRUE)


##################
# 0.72           #
##################







# Question 2 -------------------------------------------------------------------
# What is the raw (unadjusted) mean of real earnings in 1978 for treated subjects 
# minus the mean of real earnings in 1978 for untreated subjects?


lalonde %>% 
  group_by(treat) %>% 
  summarise(earnisngs = mean(re78))


earn_treat = lalonde %>% filter(treat == 1) %>% summarise(mean(re78))
earn_cont = lalonde %>% filter(treat == 0) %>% summarise(mean(re78))

earn_treat - earn_cont


##################
# -$635          #
##################








# Question 3 -------------------------------------------------------------------
# Fit a propensity score model. Use a logistic regression model, where the 
# outcome is treatment. Include the 8 confounding variables in the model as 
# predictors, with no interaction terms or non-linear terms (such as squared 
# terms). Obtain the propensity score for each subject.

# What are the minimum and maximum values of the estimated
# propensity score? 


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
             family = binomial())



#show coefficients etc
summary(model)

#create propensity score
pscore <- model$fitted.values

min(pscore)
max(pscore)


####################################
# minimum=0.009, maximum=0.85      #
####################################







# Question 4 -------------------------------------------------------------------

# Now carry out propensity score matching using the Match function. 
# 
# Before using the Match function, first do:
#   
# >set.seed(931139)
# 
# Setting the seed will ensure that you end up with a matched data set that is the 
# same as the one used to create the solutions.
# 
# Use options to specify pair matching, without replacement, no caliper. 
# 
# Match on the propensity score itself, not logit of the propensity score.  
# Obtain the standardized differences for the matched data.

# What is the standardized difference for married?


# set seed
set.seed(931139)

#do greedy matching on logit(PS) using Match without a caliper
psmatch <- Match(Tr = lalonde$treat, 
                 M = 1,
                 X = pscore,
                 caliper = NULL,
                 replace = FALSE)

matched <- lalonde[unlist(psmatch[c("index.treated","index.control")]), ]

#get standardized differences
matchedtab <- CreateTableOne(vars = xvars,
                             strata = "treat", 
                             data = matched, 
                             test = FALSE)

print(matchedtab, smd = TRUE)



###################
# 0.027           #
###################












# Question 5 -------------------------------------------------------------------
# For the propensity score matched data:
#   
# Which variable has the largest standardized difference?


# black = 0.852




###################
# black           #
###################











# Question 6 -------------------------------------------------------------------
# Re-do the matching, but use a caliper this time. Set the
# caliper=0.1 in the options in the Match function.
# 
# Again, before running the Match function, set the seed:
#   
# >set.seed(931139)
#
# How many matched pairs are there?



# set seed
set.seed(931139)

#do greedy matching on logit(PS) using Match without a caliper
psmatch <- Match(Tr = lalonde$treat, 
                 M = 1,
                 X = pscore,
                 caliper = .1,
                 replace = FALSE)

matched <- lalonde[unlist(psmatch[c("index.treated","index.control")]), ]

#get standardized differences
matchedtab <- CreateTableOne(vars = xvars,
                             strata = "treat", 
                             data = matched, 
                             test = FALSE)

print(matchedtab, smd = TRUE)



###################
# 111             #
###################








# Question 7 -------------------------------------------------------------------
# Use the matched data set (from propensity score matching with caliper=0.1) 
# to carry out the outcome analysis. 
#
# For the matched data, what is the mean of real earnings in 1978 for treated 
# subjects minus the mean of real earnings in 1978 for untreated subjects? 



matched %>% 
  group_by(treat) %>% 
  summarise(earnisngs = mean(re78))


#outcome analysis
earn_treat = matched %>% filter(treat == 1) %>% summarise(mean(re78))
earn_cont = matched %>% filter(treat == 0) %>% summarise(mean(re78))

#pairwise difference
diffy <- earn_treat - earn_cont
diffy




###################
# $1246.81        #
###################









# Question 8 -------------------------------------------------------------------
# Use the matched data set (from propensity score matching with caliper=0.1) 
# to carry out the outcome analysis.
# 
# Carry out a paired t-test for the effect of treatment on earnings. What are
# the values of the 95% confidence interval?



#paired t-test
t.test(re78 ~ treat, data = matched, paired = TRUE)


#########################
# (-420.03, 2913.64)    #
#########################






















