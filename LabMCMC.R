##Install required packages (cr. stackOverflow)
list.of.packages <- c("ggplot2", "MCMCpack", "MASS", "coda")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if (length(new.packages)) install.packages(new.packages)
require(MCMCpack)
require(MASS)
require(ggplot2)
require(coda)

#data taken from BIO 249, a course at HSPH
cnames <- c("npreg", "glu", "bp", "skin", "bmi", "ped", "age", "diabetes")
diab <- read.table("azdiabetes.txt",header = TRUE, sep = "", col.names = cnames)
data(birthwt)

##################################
###Begin Gibbs Sampler Examples###
##################################

##First we implement a Gibbs sampler using default settings
gibbsPosterior1 <- MCMCregress(bmi ~ glu + bp + age, b0 = 0, B0 = 0.1, sigma.mu = 5, sigma.var = 25, data = diab)
hist(gibbsPosterior1[,"glu"])
hist(gibbsPosterior1[,"bp"])
hist(gibbsPosterior1[,"age"])
coda::traceplot(gibbsPosterior1)
coda::acfplot(gibbsPosterior1)
coda::HPDinterval(gibbsPosterior1)
plot(gibbsPosterior1)

###Next, we set the prior for the estimates of the Betas to 10 and look at the differences
gibbsPosterior2 <- MCMCregress(bmi ~ glu + bp + age, b0 = 10, B0 = 0.1, sigma.mu = 5, sigma.var = 25, data = diab)
hist(gibbsPosterior2[,"glu"])
hist(gibbsPosterior2[,"bp"])
hist(gibbsPosterior2[,"age"])
coda::traceplot(gibbsPosterior2)
coda::acfplot(gibbsPosterior2)
coda::HPDinterval(gibbsPosterior2)


###Here, we manually set the burn in and the the amount of thinning we will do to the data.
gibbsPosterior3 <- MCMCregress(bmi ~ glu + bp + age, b0 = 0, B0 = 0.1, sigma.mu = 5, sigma.var = 25, data = diab, burnin = 1000, thin = 80)
hist(gibbsPosterior3[,"glu"])
hist(gibbsPosterior3[,"bp"])
hist(gibbsPosterior3[,"age"])
coda::traceplot(gibbsPosterior3)
coda::acfplot(gibbsPosterior3)
coda::HPDinterval(gibbsPosterior3)

df.gibbs1 <- data.frame("bp" = gibbsPosterior1[,"bp"])
df.gibbs2 <- data.frame("bp" = gibbsPosterior2[,"bp"])

df.gibbs1$Group <- "Default"
df.gibbs2$Group <- "UpdatedPrior"

bothGibbs <- rbind(df.gibbs1, df.gibbs2)
colnames(bothGibbs)[1] <- "bp"
gibbsGraph1 <- ggplot(bothGibbs,aes(x = bp,color = Group)) + geom_density() + facet_grid(Group ~ .)
gibbsGraph1

##########################################
####Begin Metropolis-Hastings Examples####
##########################################

#Default Uniform Prior
mh1 <- MCMClogit(low ~ age + as.factor(race) + smoke, data = birthwt, verbose = 1000)
plot(mh1)
summary(mh1)

mh1Summ <- summary(mh1)
mh1Summ$Prior <- "Uniform"

#Default MVN Prior
mh2 <- MCMClogit(low ~ age + as.factor(race) + smoke, b0 = 0, B0 = .001, data = birthwt)
plot(mh2)
summary(mh2)

mh2Summ <- summary(mh2)
mh2Summ$Prior <- "MVN"


###The following code is an example used by the creators of the CRAN package MCMC Pack
###This shows how we can have a user defined proprosal distribution
logpriorfun <- function(beta) {
  sum(dcauchy(beta, log = TRUE))
}

mh3 <- MCMClogit(low ~ age + as.factor(race) + smoke,
  data = birthwt, user.prior.density = logpriorfun,
  logfun = TRUE)
plot(mh3)
summary(mh3)

df.mh1 <- data.frame("smoke" = mh1[, c(1:5)])
df.mh2 <- data.frame("smoke" = mh2[, c(1:5)])
df.mh3 <- data.frame("smoke" = mh3[, c(1:5)])
df.mh1$Group <- "Uniform"
df.mh2$Group <- "MVN"
df.mh3$Group <- "Cauchy"
allMH2 <- rbind(df.mh1, df.mh2, df.mh3)

colnames(allMH2) <- c("intercept", "age", "race2", "race3", "smoke", "Group")
mhGraphSmoke <- ggplot(allMH2,aes(x = smoke ,color = Group)) + geom_density() + facet_grid(Group ~ .)
mhGraphSmoke

