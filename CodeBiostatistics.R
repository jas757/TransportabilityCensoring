load("Data.rda")
################################################################################################################
############## This code runs on simulated data that is NOT the NLST and NHANES data
############## it outputs the doubly robust estimator for the Brier risk
################################################################################################################


# Loading libraries
set.seed(1)
library(rpart)
library(survival)
library(party)
library(randomForestSRC)
library(glmnet)
#library(mboost)
library(MASS)
library(pec)
library(mgcv)

time.point = 1825/365


# Split NLST into a training and a test set
nlst = combined.data[combined.data$Sk == 1, ]
nhanes =  combined.data[combined.data$Sk == 0, ]

temp = sample(1:sum(combined.data$Sk ==1), sum(combined.data$Sk ==1)/2, replace = FALSE)
temp.not = !(1:sum(combined.data$Sk ==1) %in% temp)
nlst.train = nlst[temp, ]
nlst.test = nlst[temp.not, ]

combined.test = rbind(nlst.test, nhanes)

##########################################################################################
########## Fit a random forest model
##########################################################################################

rand.for.mod = rfsrc(Surv(obs, delta)~.-Sk-weights, data = nlst.train)

# Getting the Survival Curves. 
pred.rf <- predict(rand.for.mod, newdata = nlst.test, proximity = FALSE)

# predsurvAFT compute P(T>t|W) tval is t, fit is the AFT fit and zval is the covariate values.
predsurvRF = function(time.point, cov.index){
time.used = sum(pred.rf$time.interest < time.point) + 1
surv = c(1, pred.rf$survival[cov.index, ])[time.used]
return(surv)
}

pred.x = rep(NA, nrow(nlst.test))
for(i in 1:nrow(nlst.test)){
pred.x[i] = predsurvRF(time.point, i)
}

##########################################################################################
########## Evaluate the model
##########################################################################################

# Create the data
obs = combined.test$obs
delta = combined.test$delta
xx = combined.test[ ,3:16]
S = combined.test$Sk
weights.samp = combined.test$weights

# Estimator of the censoring probabilities using a certain form of truncation.
GfuncKM=function(obs,delta,dtype)
{ 
  n <- length(obs)
  # Changing the dataset to account for truncation. 
  aa=datach(obs,delta,dtype)
  # Observed time after truncation
  obs=aa[[1]]
  # Failure indicator after truncation.
  delta=aa[[2]]
  
  #Calculating the hazard estimator. 
  hazC=mapply(function(xx,dd){dd/sum((obs>=xx))},xx=obs,dd=1-delta)
  surC_km=mapply(function(xx){prod(1-hazC[obs<=xx])},xx=obs)
  return(list(surC_km,obs,delta))
}



# Calculates The estimator for the censoring distribution using a survival tree
# Input: obs = observed time, delta = failure indicator, dtype = truncation level and method.
# Output: List of bar -(tilde T|W), failure indicator changed to account for truncation
# and observed time also changed to account for truncation.
GfuncSurvivalTree=function(obs,delta,dtype,xx)
{

num = length(obs)
nu = num

# Fitting the Cox model
# Creating the data frame
data.used <- data.frame(obs, 1 - delta, xx)
names(data.used)[1:2] = c("obs", "delta.cens")
surv.tree = rpart(Surv(obs,delta.cens)~., data = data.used, minbucket = 30, maxdepth = 1)

# Getting the Survival Curves. 
pred.surv.tree <- predict(surv.tree, proximity = FALSE)

# Finding the terminal nodes
sett=unique(surv.tree[['where']])
nset=length(sett)

cens.est = matrix(0, nrow = num, ncol = num)
obs.used = rep(NA, num)
delta.used = rep(NA, num)

for (i in 1:nset){
# Finding the subset corresponding the ith node of the tree
subset=(1:nu)[surv.tree[['where']]==sett[i]]
nlen=length(subset)
# Observed time within the node
sobs=obs[subset]
# Failure indicators within each node.
sdelta=delta[subset]

# Doing truncation within that subset
# Changing the dataset to account for truncation. 
aa=datach(sobs,sdelta,dtype = dtype)
# Observed time after truncation
sobs=aa[[1]]
# Failure indicator after truncation.
sdelta = aa[[2]]

obs.used[subset] = sobs
delta.used[subset] = sdelta

# Calculating the KM estimator censoring curve within a node
# Calculating the jumps in the KM estimator
hazC=mapply(function(xx,dd){dd/sum((sobs>=xx))},xx=sobs,dd=1-sdelta)
surC_km=mapply(function(xx){prod(1-hazC[sobs<=xx])},xx=obs)
cens.est[subset, ] = matrix(surC_km,nrow=length(subset),ncol=length(surC_km),byrow=TRUE)
}

return(list(cens.est, obs.used, delta.used, surv.tree[['where']]))
}



# Changes the dataset to account for truncation
# Input: obs = observed time, delta = failure indicator, dtype = how truncation is done.
# b is method 2 and a is method 1 and the number afterwards indicates what the truncation level is.
# Output: obs = observed failure time adjusted for truncation, delta = failure indicator adjusted for truncaiton.
datach=function(obs,delta,dtype)
{
  nu = length(obs)
  if(dtype=="b0")
  {
    delta[obs==max(obs)]=TRUE
  }
  
  
  
  if(dtype=="b5")
  {
    delta[order(obs)][floor(nu-0.05*nu):nu]=TRUE
  }
  
  
  if(dtype=="b10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=TRUE
  }
  
  
  if(dtype=="b15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=TRUE
  }
  
  if(dtype=="a5")
  {    
    delta[order(obs)][floor(nu-0.05*nu):nu]=T
    obs[order(obs)][floor(nu-0.05*nu):nu]=obs[order(obs)][floor(nu-0.05*nu)]    
  } 
  if(dtype=="a10")
  {
    delta[order(obs)][floor(nu-0.10*nu):nu]=T
    obs[order(obs)][floor(nu-0.10*nu):nu]=obs[order(obs)][floor(nu-0.10*nu)]    
  } 
  
  
  if(dtype=="a15")
  {
    delta[order(obs)][floor(nu-0.15*nu):nu]=T
    obs[order(obs)][floor(nu-0.15*nu):nu]=obs[order(obs)][floor(nu-0.15*nu)]    
  }   
  
  if(dtype=="none")
  {
    obs=obs
    delta=delta
  }   
  return(list(obs,delta)) 
} 




# Calculating the conditional expectation
# Input: obs = observed time, delta = failure indicator, xx = vector of covariates
# mtype = what type of cond exp is being calculated. 
# Output: Matrix of conditional expectations m_{1i}(tilde T_j) for j, i = 1, ldots ,n
# time point is the t in the equation m1[j, i] = P(T >t|T >T_j, W_i)
mfuncBrier=function(obs,delta,xx, mtype, time.point)
{
  num=length(obs)
  
  if(mtype=="cox")
  {
    m1=coxxBrier(obs,delta, xx, time.point)
  }

  if(mtype=="weibull")
  {
    m1=AftWeibullBrier(obs,delta,xx, time.point)
  }

  if(mtype=="lognorm")
  {
    m1=AftLogNormBrier(obs,delta,xx, time.point)
  }

  if(mtype=="loglog")
  {
    m1 = AftLogLogBrier(obs,delta,xx, time.point)
  }

  if(mtype=="rand.for")
  {
    m1 = randomForestBrier(obs,delta,xx, time.point)
  }
 
  return(m1)  
}



####################################################################################################
###### The following functions calcualate the estimator for P(T >t|W)
#####################################################################################################

# Calculating the model for the Conditional Expectations using the Cox Model.
coxxBrier = function(obs,delta,xx, time.point){

# Number of observations
n = length(delta)

# Defininf the matrix where m1[j, i] = P(T>t|T >T_j, W_i)
m1 <- matrix(1, ncol = n, nrow = n)

data.nn = data.frame(obs, delta, xx)
# Fitting the censoring cox model
cox.model <- coxph(Surv(obs, delta) ~ xx, data = data.nn, x = TRUE)

# Calculating the cox model
#m1[j, i] <- P(T >t|T > T_j, W_i)
for(j in 1:n){
if(delta[j] == 0 & obs[j] < time.point){
# Calculate the denominator and the numerator in the desired probability
prob.est.num = predictSurvProb(cox.model, newdata = data.nn, times = time.point)
prob.est.den = predictSurvProb(cox.model, newdata = data.frame(xx), times = obs[j])
m1[j, ] = prob.est.num/prob.est.den
}
}

# Return the probabilities
return(m1)
}


# Calculating the model for the Conditional Expectations using the random forest algorithm.
randomForestBrier = function(obs,delta,xx, time.point){

n = length(obs)
nu = length(obs)

# Fitting the Cox model
# Creating the data frame
data.used <- data.frame(obs, delta, xx)
rand.for = rfsrc(Surv(obs,delta)~ ., data = data.used, importance = "none")
  
# Getting the Survival Curves. 
pred.rf <- predict(rand.for, proximity = FALSE)

# Calculating the cox model
#m1[j, i] <- P(T >t|T > T_j, W_i)
#P(T > tau|T >u,W) = P(T>tau|W)/P(T>u|W)

# predsurvAFT compute P(T>t|W) tval is t, fit is the AFT fit and zval is the covariate values.
predsurvRF = function(time, cov.index){
time.point = sum(pred.rf$time.interest < time) + 1
surv = c(1, pred.rf$survival[cov.index, ])[time.point]
return(surv)
}

# Defininf the matrix where m1[j, i] = P(T>t|T >T_j, W_i)
m1 <- matrix(0, ncol = n, nrow = n)

# Calculating the cox model
#m1[j, i] <- P(T >t|T > T_j, W_i)
for(j in 1:n){
if(delta[j] == 0 & obs[j] < time.point){
for(i in 1:n){
# Calculate the denominator and the numerator in the desired probability
prob.est.den = predsurvRF(obs[j], i)
prob.est.num = predsurvRF(time.point, i)
if(prob.est.den != 0){
m1[j, i] = prob.est.num/prob.est.den
}
if(prob.est.den == 0){
  m1[j, i] = 0.5
}

}
}
}

# Return the probabilities
return(m1)
}


### Calculate an AFT model

# Calculating the model for the Conditional Expectations using the AFT model with a Weibull distribution.
AftWeibullBrier = function(obs,delta, xx, time.point){
    
# Number of observations
n = length(delta)

data.fit = data.frame(obs,delta, xx)
weibullfit = survreg(Surv(obs, delta) ~ ., dist = "weibull", data = data.fit, maxiter = 1)

# predsurvAFT compute P(T>t|W) tval is t, fit is the AFT fit and zval is the covariate values.
predsurvAFT=function(tval,zval,fit){
z=as.vector(zval)
lp=fit[['coefficients']]%*% unlist(as.vector(c(1,z))) # alpha’z
ztheta=(log(tval)-lp)/fit[['scale']] # z(theta) from Slide 1
surv=exp(-exp(ztheta))
return(surv)
}

# Defininf the matrix where m1[j, i] = P(T>t|T >T_j, W_i)
m1 <- matrix(1, ncol = n, nrow = n)

# Calculating the predictions from AFT
#m1[j, i] <- P(T >t|T > T_j, W_i)
for(j in 1:n){
if(delta[j] == 0 & obs[j] < time.point){
for(i in 1:n){
# Calculate the denominator and the numerator in the desired probability
prob.est.den = predsurvAFT(obs[j], xx[i,], weibullfit)
prob.est.num = predsurvAFT(time.point, xx[i,], weibullfit)
m1[j, i] = prob.est.num/prob.est.den
}
}
}

return(m1)
}



# Calculating the model for the Conditional Expectations using log normal distribution.
AftLogNormBrier = function(obs,delta,xx, time.point){
    
# Number of observations
n = length(delta)

data.fit = data.frame(obs,delta, xx)
loglogfit = survreg(Surv(obs, delta) ~ xx, dist = "lognormal", data = data.fit)

# predsurvAFT compute P(T>t|W) tval is t, fit is the AFT fit and zval is the covariate values.
predsurvAFT=function(tval,zval,fit){
z=as.vector(zval)
lp=fit[['coefficients']]%*% as.vector(c(1,z)) # alpha’z
ztheta=(log(tval)-lp)/fit[['scale']] # z(theta) from Slide 1
surv=1-pnorm(ztheta)
return(surv)
}

# Defininf the matrix where m1[j, i] = P(T>t|T >T_j, W_i)
m1 <- matrix(1, ncol = n, nrow = n)

# Calculating the cox model
#m1[j, i] <- P(T >t|T > T_j, W_i)
for(j in 1:n){
if(delta[j] == 0 & obs[j] < time.point){
for(i in 1:n){
# Calculate the denominator and the numerator in the desired probability
prob.est.den = predsurvAFT(obs[j], xx[i,], loglogfit)
prob.est.num = predsurvAFT(time.point, xx[i,], loglogfit)
m1[j, i] = prob.est.num/prob.est.den
}
}
}

return(m1)
}


# Calculating the model for the Conditional Expectations using the a log-log distribution
AftLogLogBrier = function(obs,delta,xx, time.point){
    
# Number of observations
n = length(delta)
data.fit = data.frame(obs,delta, xx)
loglogfit = survreg(Surv(obs, delta) ~ xx, dist = "loglogistic", data = data.fit)

# predsurvAFT compute P(T>t|W) tval is t, fit is the AFT fit and zval is the covariate values.
predsurvAFT=function(tval,zval,fit){
z=as.vector(zval)
lp=fit[['coefficients']]%*% as.vector(c(1,z)) # alpha’z
ztheta=(log(tval)-lp)/fit[['scale']] # z(theta) from Slide 1
surv=1/(1+exp(ztheta))
return(surv)
}

# Defininf the matrix where m1[j, i] = P(T>t|T >T_j, W_i)
m1 <- matrix(1, ncol = n, nrow = n)

# Calculating the cox model
#m1[j, i] <- P(T >t|T > T_j, W_i)
for(j in 1:n){
if(delta[j] == 0 & obs[j] < time.point){
for(i in 1:n){
# Calculate the denominator and the numerator in the desired probability
prob.est.den = predsurvAFT(obs[j], xx[i,], loglogfit)
prob.est.num = predsurvAFT(time.point, xx[i,], loglogfit)
m1[j, i] = prob.est.num/prob.est.den
}
}
}

return(m1)
}



###########################################################################################
####### externalBrierTreeKM creates imputation using survival tree for G_0 and mtype 
####### determines the model for survival curve
###########################################################################################

# external calculates the parameter vector using the Brier loss using a tree model
# Inputs: obs = observed time, delta = failure indicator, x1-x5 covariates,
# mtype = type of conditional expectation, dtype = truncation method and level.
externalBrierTreeKM = function(obs,delta,xx, mtype,dtype, time.point, pred.x)
{
  n = length(obs)
  nu = length(obs)

  # Creating the new T(t) dataset
  obs.t = pmin(obs, time.point)
  delta.t = delta * (obs <= time.point) + (obs > time.point)
  n = length(obs)
  
  # Calculating the conditional expectation    
  m1=mfuncBrier(obs,delta,xx, mtype, time.point)

  # Calculating the conditional censoring distribution.
  tem = GfuncSurvivalTree(obs.t,delta.t,dtype, xx)
  # Calculating the censoring distribution
  surC_rf=tem[[1]]
  # Finding which observations fall in which terminal node
  term.node = tem[[4]]

  # Calculating a0, a1, b0, b1, c0, c1
  a0=delta.t/diag(surC_rf)
  a1 = a0 * ((obs > time.point) - pred.x)^2
  
  b0=(1-delta.t)/diag(surC_rf)
  b1=b0 * (diag(m1) - 2 * pred.x * diag(m1) + pred.x^2)

c0 <- rep(NA, n)
c1 <- rep(NA, n)

# Finding the terminal nodes
sett=unique(term.node)
nset=length(sett)

for (i in 1:nset){
# Finding the subset corresponding the ith node of the tree
subset=(1:nu)[term.node==sett[i]]
nlen=length(subset)
# Observed time within the node
sobs=obs.t[subset]
# Failure indicators within each node.
sdelta=delta.t[subset]


  kk=mapply(function(tt){sum((tt<=sobs))},tt=sobs)
  c0[subset]=mapply(function(tt){sum(b0[subset]*(sobs<=tt)/kk)},tt=sobs)
  c1[subset]=mapply(function(tt,i){sum(b0[subset]*(sobs<=tt)*(m1[subset,i] - 2 * pred.x[subset] * m1[subset,i] + pred.x[subset]^2)/kk)},tt=sobs,i=1:nlen)
}

dr.imp = a1 + b1 - c1
bj.imp = delta.t * (as.numeric(obs > time.point) - pred.x)^2 + (1 - delta.t) * (diag(m1) - 2 * pred.x * diag(m1) + pred.x^2)
ipw.imp = a1

parms = c(a0,a1,b0,b1,c0,c1,obs,delta,diag(m1),nu)
  
  return(list(dr.imp = dr.imp, bj.imp = bj.imp, ipw.imp = ipw.imp))
}



#########################################################################################################
#############  Create the doubly robust transportability estimator using an AFT model with a 
############   Weibull distribution
#########################################################################################################

mtype = "rand.for"
dtype = "b10"
# Create survival imputations 
imp.all = externalBrierTreeKM(obs[S==1],delta[S==1], xx[S==1, ], mtype, dtype, time.point, pred.x[S==1]) 
imp.dr = imp.all$dr.imp

# Regress the pseudo outcomes on X
data.lm = data.frame(imp.dr, xx[S==1, ])
lin.mod = gam(imp.dr ~ s(age) + s(bmi) + race + educat2 + smokeyr + gender + marital + pkyr + smokelive + diagadas +diagdiab + diagemph+ diaghear + diaghype, data = data.lm)
h.x = pmax(predict(lin.mod, newdata = data.frame(xx)), 0)

# Create source population model
data.source.mod = data.frame(S, xx, weights.samp)
s.mod = glm(S~.-weights.samp, data = data.source.mod, family = "binomial", weights = weights.samp)
p.x = predict(s.mod, newdata = data.frame(xx), type = "response")

# Create alpha variable for dr estimator
alpha.var = rep(NA, nrow(xx))
alpha.var[S==0] = weights.samp[S==0] * h.x[S==0]/sum(weights.samp[S==0]) 
alpha.var[S==1] = (1 - p.x[S==1])/p.x[S==1] * (imp.dr - h.x[S==1])/sum(weights.samp[S==0])
dr.brier.targ = sum(alpha.var)
