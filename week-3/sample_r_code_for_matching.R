###################
#RHC Example

#install packages
# install.packages("tableone")
# install.packages("Matching")

# install.packages("optmatch")

#load packages
library(tableone)
library(Matching)
library(MatchIt)


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

#outcome analysis
y_trt <- mydata$died[mydata$treatment==1]
y_con <- mydata$died[mydata$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)



############################################
#do greedy matching on Mahalanobis distance
############################################

greedymatch <- Matching::Match(Tr = treatment,
                               M = 1,
                               X = mydata[xvars],
                               replace = FALSE)

matched <- mydata[unlist(greedymatch[c("index.treated","index.control")]), ]


#get table 1 for matched data with standardized differences
matchedtab0 <- tableone::CreateTableOne(vars = xvars, 
                                        strata ="treatment",
                                        data = matched, 
                                        test = T)

print(matchedtab0, smd = TRUE)

#outcome analysis
y_trt <- matched$died[matched$treatment==1]
y_con <- matched$died[matched$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)

#McNemar test
table(y_trt,y_con)

mcnemar.test(matrix(c(973,513,395,303),2,2))



############################################
# other methods to match - mahalanobis
############################################


match1 <- MatchIt::matchit(treatment ~ 
                             ARF + CHF + Cirr + colcan + Coma + lungcan +
                             MOSF + sepsis + age + female + meanbp1, 
                           data = mydata,
                           replace = TRUE,
                           distance = "mahalanobis")

summary(match1)
plot(summary(match1))
dfm1 <- MatchIt::match.data(match1)

#get table 1 for matched data with standardized differences
matchedtab1 <- tableone::CreateTableOne(vars = xvars, 
                                        strata ="treatment",
                                        data = dfm1, 
                                        test = T)

print(matchedtab1, smd = TRUE)


#outcome analysis
y_trt <- dfm1$died[dfm1$treatment==1]
y_con <- dfm1$died[dfm1$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)



############################################
# other methods to match - nearest glm
############################################


match2 <- MatchIt::matchit(treatment ~ 
                             ARF + CHF + Cirr + colcan + Coma + lungcan +
                             MOSF + sepsis + age + female + meanbp1, 
                           data = mydata,
                           replace = TRUE,
                           method = "nearest",
                           distance = "glm")

summary(match2)
plot(summary(match2))
dfm2 <- MatchIt::match.data(match2)

#get table 1 for matched data with standardized differences
matchedtab2 <- tableone::CreateTableOne(vars = xvars, 
                                        strata ="treatment",
                                        data = dfm2, 
                                        test = T)

print(matchedtab2, smd = TRUE)


#outcome analysis
y_trt <- dfm2$died[dfm2$treatment==1]
y_con <- dfm2$died[dfm2$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)




############################################
# other methods to match - full glm probit
############################################


match3 <- MatchIt::matchit(treatment ~ 
                             ARF + CHF + Cirr + colcan + Coma + lungcan +
                             MOSF + sepsis + age + female + meanbp1, 
                           data = mydata,
                           method = "full",
                           distance = "glm",
                           link = "probit")

summary(match3)
plot(summary(match3))
dfm3 <- MatchIt::match.data(match3)

#get table 1 for matched data with standardized differences
matchedtab3 <- tableone::CreateTableOne(vars = xvars, 
                                        strata ="treatment",
                                        data = dfm3, 
                                        test = T)

print(matchedtab3, smd = TRUE)


#outcome analysis
y_trt <- dfm3$died[dfm3$treatment==1]
y_con <- dfm3$died[dfm3$treatment==0]

#pairwise difference
diffy<-y_trt-y_con

#paired t-test
t.test(diffy)





plot(match3, type = "jitter")
plot(match3, type = "hist")






##########################
#propensity score matching
#########################

mydata

#fit a propensity score model. logistic regression

psmodel<-glm(treatment~ARF+CHF+Cirr+colcan+Coma+lungcan+MOSF+
               sepsis+age+female+meanbp1,
             family=binomial(),data=mydata)

#show coefficients etc
summary(psmodel)
#create propensity score
pscore<-psmodel$fitted.values


#do greedy matching on logit(PS) using Match without a caliper
logit <- function(p) {log(p)-log(1-p)}
psmatch<-Match(Tr=mydata$treatment,M=1,X=logit(pscore),replace=FALSE)
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





# use matchit for propensity score, nearest neighbor matching
m.out <- MatchIt::matchit(treatment ~
                            ARF + CHF + Cirr + colcan + Coma + lungcan +
                            MOSF + sepsis + age + female + meanbp1,
                          data = mydata,
                          method = "nearest")

summary(m.out)
plot(m.out, type = "jitter")
plot(m.out, type = "hist")

































