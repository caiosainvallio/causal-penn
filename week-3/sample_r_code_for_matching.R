###################
#RHC Example

#install packages
install.packages("tableone")
install.packages("Matching")

#load packages
library(tableone)
library(Matching)


# read in data
load(url("https://biostat.app.vumc.org/wiki/pub/Main/DataSets/rhc.sav"))
#view data
View(rhc)

rhc <- tibble::tibble(rhc)

#treatment variables is swang1
#x variables that we will use
#cat1: primary disease category
#age
#sex
#meanbp1: mean blood pressure

#create a data set with just these variables, for simplicity
ARF <- as.numeric(rhc$cat1 == 'ARF')
CHF <- as.numeric(rhc$cat1 == 'CHF')
Cirr <- as.numeric(rhc$cat1 == 'Cirrhosis')
colcan <- as.numeric(rhc$cat1 == 'Colon Cancer')
Coma <- as.numeric(rhc$cat1 == 'Coma')
COPD <- as.numeric(rhc$cat1 == 'COPD')
lungcan <- as.numeric(rhc$cat1 == 'Lung Cancer')
MOSF <- as.numeric(rhc$cat1 == 'MOSF w/Malignancy')
sepsis <- as.numeric(rhc$cat1 == 'MOSF w/Sepsis')
female <- as.numeric(rhc$sex == 'Female')
died <- as.numeric(rhc$death == 'Yes')
age <- as.numeric(rhc$age)
treatment <- as.numeric(rhc$swang1 == 'RHC')
meanbp1 <- as.numeric(rhc$meanbp1)

#new dataset
mydata <- tibble::tibble(ARF,CHF,Cirr,colcan,Coma,lungcan,MOSF,sepsis,
                       age,female,meanbp1,treatment,died)

#covariates we will use (shorter list than you would use in practice)
xvars <- c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#look at a table 1
table1 <- tableone::CreateTableOne(vars = xvars,
                                  strata = "treatment", 
                                  data = mydata,
                                  test = T)
## include standardized mean difference (SMD)
print(table1,smd=TRUE)






############################################
#do greedy matching on Mahalanobis distance
############################################

greedymatch <- Matching::Match(Tr = treatment,
                               M = 1,
                               X = mydata[xvars],
                               replace = FALSE)

matched <- mydata[unlist(greedymatch[c("index.treated","index.control")]), ]


#get table 1 for matched data with standardized differences
matchedtab1 <- tableone::CreateTableOne(vars = xvars, 
                                        strata ="treatment",
                                        data = matched, 
                                        test = T)

print(matchedtab1, smd = TRUE)



#outcome analysis
y_trt <- matched$died[matched$treatment==1]
y_con <- matched$died[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)
t.test(died ~ treatment)

# Point estimate: 0.0540293
# Difference in probability of death if everyone received  RHC versus
# no one received RHC is 0.05 (i.e., higher risk of death in RHC group)
# 95%  CI: (0.027, 0.080)
# P-value: <0.001



#McNemar test
table(y_trt,y_con)

mcnemar.test(matrix(c(973,513,395,303),2,2))



##########################
#propensity score matching
#########################

#fit a propensity score model. logistic regression

psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1+aps,
             family=binomial(),data=mydata)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values


#do greedy matching on logit(PS) using Match with a caliper

logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE,caliper=.2)
matched<-mydata[unlist(psmatch[c("index.treated","index.control")]), ]
xvars<-c("ARF","CHF","Cirr","colcan","Coma","lungcan","MOSF","sepsis",
         "age","female","meanbp1")

#get standardized differences
matchedtab1<-CreateTableOne(vars=xvars, strata ="treatment", 
                            data=matched, test = FALSE)
print(matchedtab1, smd = TRUE)

#outcome analysis
y_trt<-matched$died[matched$treatment==1]
y_con<-matched$died[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)
