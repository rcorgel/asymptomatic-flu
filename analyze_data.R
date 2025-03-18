################################################################################
# File Name: analyze_data                                                      #
#                                                                              #
# Purpose:   Analyze asymptomatic influenza log percents with Bayesian linear  #
#            regression.                                                       #
# Steps:                                                                       # 
#            1. Set-up script                                                  #
#            2. Load and clean data                                            #
#            3. Transform variables                                            #
#            4. Perform regression analysis                                    #
#                                                                              #
# Project:   Asymptomatic Influenza Associations                               #
################################################################################

####################
# 1. SET-UP SCRIPT #
####################

# Start with a clear environment
rm(list = ls())

# Load libraries
library(tidyverse)
library(lubridate)
library(assertr)
library(mvnfast)
library(corrplot)
library(MCMCpack)
library(car)
library(Hmisc)
library(cowplot)
library(reshape2)
library(moments)

# Set the seed
set.seed(12345)

# Set the directory
setwd('/Users/rcorgel/Documents/asymptomatic-flu-proj/')

##########################
# 2. LOAD AND CLEAN DATA #
##########################

# Load combined data
data <- read.csv('tmp/analysis_dat.csv')

# Keep complete cases
analysis_dat <- data[complete.cases(data), ]
# 81.2% of data remains after dropping
length(analysis_dat$county_fips) / length(data$county_fips)
# Overall 2078 out of 3144 counties have complete and un-censored data (66.0%)

# Remove index variable
analysis_dat <- analysis_dat |> dplyr::select(-(X))

##########################
# 3. TRANSFORM VARIABLES #
##########################

# Log transform asymptomatic flu for regression
analysis_dat$percent_asymp_flu_log <- log(analysis_dat$percent_asymp_flu * 100)
analysis_dat$percent_asymp_flu_scale <- analysis_dat$percent_asymp_flu * 100

# Change disease counts to percents (over all patients in season)
analysis_dat$percent_flu <- analysis_dat$total_flu / analysis_dat$total_all_cause * 1000
analysis_dat$percent_rsv <- analysis_dat$total_rsv / analysis_dat$total_all_cause * 1000
analysis_dat$percent_covid <- analysis_dat$total_covid / analysis_dat$total_all_cause * 1000

# Healthcare access per 100,000 population
analysis_dat$rate_icu_beds <- analysis_dat$icu_beds / analysis_dat$population * 100000
analysis_dat$rate_doctors <- analysis_dat$rate_doctors * 100000 # already a percent

# Scale down median hh income and population
analysis_dat$median_hh_income_scale <- analysis_dat$median_hh_income / 10000 # income in 10,000
analysis_dat$population_scale <- analysis_dat$population / 100000 # population in 100,000

# Scale all remaining percent variables to be 0-100 for interpretation
analysis_dat$percent_65_plus_flu_scale <- analysis_dat$percent_65_plus_flu * 100
analysis_dat$percent_white_flu_scale <- analysis_dat$percent_white_flu * 100
analysis_dat$percent_female_flu_scale <- analysis_dat$percent_female_flu * 100
analysis_dat$percent_poor_health_scale <- analysis_dat$percent_poor_health * 100
analysis_dat$percent_obese_scale <- analysis_dat$percent_obese * 100
analysis_dat$percent_mental_distress_scale <- analysis_dat$percent_mental_distress * 100
analysis_dat$percent_rural_scale <- analysis_dat$percent_rural * 100
analysis_dat$percent_vaccinated_scale <- analysis_dat$percent_vaccinated * 100

##################################
# 4. PERFORM REGRESSION ANALYSIS #
##################################

# Limit data to regression variables
regress_dat <- analysis_dat[, c(1, 19, 20, 24, 33, 34, 35, 36, 37, 38, 39,
                 41, 42, 43, 45, 46, 47, 48)]

# Save regression data
write.csv(analysis_dat, "out/regress_dat.csv")
write.csv(analysis_dat, "asymptomatic-flu/data/regress_dat.csv")

# Examine variable correlations
corr <- cor(regress_dat[, c(-1)])
corrplot(corr, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

# Quick frequentist testing
model <- lm(percent_asymp_flu_scale ~ percent_65_plus_flu_scale + 
              percent_female_flu_scale + percent_vaccinated_scale +
              percent_obese_scale + percent_mental_distress_scale + 
              percent_white_flu_scale + 
              avg_humidity + avg_temp + 
              rate_doctors + rate_icu_beds + median_hh_income_scale +
              percent_flu + percent_rsv + percent_covid +
              percent_rural_scale, 
            data = regress_dat)
summary(model)
AIC(model)
vif(model) # quick multicollinearity check
plot(model) # quick check of linear reg assumptions

########################################
# Variable selection: All combinations #
########################################

# Create a list of all possible variable combinations (16)
combo_list <- NULL
for (i in 1:15) {
  combos <- as.matrix(t(combn(2:16, i)))
  combos_full <- cbind(combos, matrix(data = NA, nrow = nrow(combos), ncol = 15-i))
  combo_list <- rbind(combo_list, combos_full)
}
nrow(combo_list) # over 32,000 variable combinations

# Set the number of samples
B <- 10000 
# Define the number of observations
n <- nrow(regress_dat)
# Define the outcome vector
Y <- regress_dat$percent_asymp_flu_scale

# Specify the design matrix
X <- model.matrix(percent_asymp_flu_scale ~ 
                    # immunity variables
                    percent_65_plus_flu_scale + percent_female_flu_scale + 
                    percent_vaccinated_scale + percent_obese_scale + 
                    percent_mental_distress_scale + 
                    # biological variables
                    avg_humidity + avg_temp + 
                    # healthcare access variables
                    rate_doctors + rate_icu_beds + 
                    # SES variables
                    median_hh_income_scale + percent_white_flu_scale +
                    # healthcare burden variables
                    percent_flu + percent_rsv + percent_covid +
                    # urbanicity variable
                    percent_rural_scale, 
                  data = regress_dat)

# Set an empty vector for WAIC values
all_waic <- vector(length = 35000)

# Loop through all variable combinations, calculating WAIC
for (i in 1:nrow(combo_list)) {
  # Print a tracking variable
  print(i)
  
  # Re-specify the design matrix based on the variable combination
  X_new <- X[, c(1,c(combo_list[i,])[!is.na(c(combo_list[i,]))])]
  
  # Define the number of independent variables (new design matrix)
  K <- ncol(X_new)
  
  # Calculate regression inputs
  bhat <- c(solve(t(X_new)%*%X_new)%*%(t(X_new)%*%Y))
  SSY	<- t(Y - X_new%*%bhat)%*%(Y - X_new%*%bhat)
  XtXi <- solve(t(X_new)%*%X_new)
  
  # Create a matric to store betas
  rbeta <- matrix(0, nrow = B, ncol = K)
  
  # Sample sigma values
  rsig <- rinvgamma(B, (n-K)/2, (1/2)*SSY)
  
  # Calculate a Var/Cov matrix from each sigma and sample betas
  for (j in 1:B){
    CovX <- rsig[j]*XtXi
    rbeta[j,]	<- c(rmvn(1, mu = bhat, sigma = CovX))
  }

  # Set empty rectors for point-wise WAIC calculations
  lppdv <- vector(length = n)
  pwaicv <- vector(length = n)
  
  # Loop through all counties
  for(h in 1:n){
    Yi <- Y[h]
    mu <- rbeta %*% X_new[h, ]
    liki <- dnorm(Yi, mu, sqrt(rsig))
    lppdv[h] <- 1/B*sum(liki)
    pwaicv[h]	<- (1/(B-1))*sum((log(liki) - mean(log(liki)))^2)
  }
  
  # Calculate WAIC
  lppd	<- sum(log(lppdv))
  pwaic	<- sum(pwaicv)
  WAIC	<- -2*lppd +2*pwaic
  
  # Add WAIC value to vector
  all_waic[i] <- WAIC
  
}

# Merge on WAIC values to variable combinations
combo_list$WAIC <- all_waic

# Identify the model with the lowest WAIC
min(combo_list$WAIC)

# Re-specify the design matrix based on the variable combination lowest WAIC
X_new <- X[, c(1,c(combo_list[29527,])[!is.na(c(combo_list[29527,]))])]

# Define the number of independent variables (new design matrix)
K <- ncol(X_new)

# Calculate regression inputs
bhat <- c(solve(t(X_new)%*%X_new)%*%(t(X_new)%*%Y))
SSY	<- t(Y - X_new%*%bhat)%*%(Y - X_new%*%bhat)
XtXi <- solve(t(X_new)%*%X_new)

# Create a matric to store betas
rbeta <- matrix(0, nrow = B, ncol = K)

# Sample sigma values
rsig <- rinvgamma(B, (n-K)/2, (1/2)*SSY)

# Calculate a Var/Cov matrix from each sigma and sample betas
for (j in 1:B){
  CovX <- rsig[j]*XtXi
  rbeta[j,]	<- c(rmvn(1, mu = bhat, sigma = CovX))
}

# Set empty rectors for point-wise WAIC calculations
lppdv <- vector(length = n)
pwaicv <- vector(length = n)

# Loop through all counties
for(h in 1:n){
  Yi <- Y[h]
  mu <- rbeta %*% X_new[h, ]
  liki <- dnorm(Yi, mu, sqrt(rsig))
  lppdv[h] <- 1/B*sum(liki)
  pwaicv[h]	<- (1/(B-1))*sum((log(liki) - mean(log(liki)))^2)
}

# Calculate WAIC
lppd	<- sum(log(lppdv))
pwaic	<- sum(pwaicv)
WAIC	<- -2*lppd +2*pwaic

# Display results
rbMat <- apply(rbeta, 2, quantile, probs = c(0.5, 0.025, 0.975))
postp <- apply(rbeta > 0, 2, mean)
outMat <- cbind(t(rbMat), postp)
rownames(outMat) <- colnames(X_new)
colnames(outMat) <- c(rownames(rbMat), 'P(b > 0)')
data_frame_all <- as.data.frame(outMat)
data_frame_all$Variables <- c('Intercept' , 'Influenza Percent 65+', 
                              'Percent Vaccinated', 'Percent Mental Distress',
                              'Avg. Temperature', 'Rate Doctors',
                              'Rate ICU Beds', 'Median HH Income',
                              'Influenza Percent White', 'Rate Influenza',
                              'Percent Rural')

# Save results for plotting
write.csv(data_frame_all, file = "tmp/data_frame_all.csv")
write.csv(as.data.frame(combo_list), file = "tmp/combo_list.csv")

# Calculate model fit statistics
yrep <- rbeta %*% t(X_new) # predicted y values for 10,000 samples
Tyrep	<- matrix(0, nrow = B, ncol = 4) # empty matrix to fill in values

# Loop through samples and calculate statistics
for(b in 1:B){
  Tyrep[b,1]	<- mean(yrep[b,])
  Tyrep[b,2]	<- var(yrep[b,])
  Tyrep[b,3]	<- skewness(yrep[b,])
  Tyrep[b,4]	<- kurtosis(yrep[b,])
}
Ty1	<- mean(Y)
Ty2 <- var(Y)
Ty3	<- skewness(Y)
Ty4	<- kurtosis(Y)

# Examine how predicted values relate to observed
mean(Tyrep[,1] > Ty1)
mean(Tyrep[,2] > Ty2)
mean(Tyrep[,3] > Ty3)
mean(Tyrep[,4] > Ty4)

# Reshape data to get distribution traces and save
traces <- as.data.frame(t(yrep))
long <- reshape2::melt(traces, id.vars = c())
write.csv(long, file = "tmp/traces_dist.csv")

############################################
# Variable selection: Background knowledge #
############################################

# Immunological & Biological Explanation
# Reset B, n, Y, X, K
# Set the number of samples
B <- 10000 
# Define the number of observations
n <- nrow(regress_dat)
# Define the outcome vector
Y <- regress_dat$percent_asymp_flu_scale

# Specify the design matrix
X <- model.matrix(percent_asymp_flu_scale ~ 
                    # immunity variables
                    percent_65_plus_flu_scale + percent_female_flu_scale + 
                    percent_vaccinated_scale + percent_obese_scale + 
                    percent_mental_distress_scale + 
                    # biological variables
                    avg_humidity + avg_temp, 
                  data = regress_dat)

# Define the number of independent variables (new design matrix)
K <- ncol(X)

# Calculate regression inputs
bhat <- c(solve(t(X)%*%X)%*%(t(X)%*%Y))
SSY	<- t(Y - X%*%bhat)%*%(Y - X%*%bhat)
XtXi <- solve(t(X)%*%X)
  
# Create a matric to store betas
rbeta <- matrix(0, nrow = B, ncol = K)
  
# Sample sigma values
rsig <- rinvgamma(B, (n-K)/2, (1/2)*SSY)
  
# Calculate a Var/Cov matrix from each sigma and sample betas
for (i in 1:B){
  CovX <- rsig[i]*XtXi
  rbeta[i,]	<- c(rmvn(1, mu = bhat, sigma = CovX))
}
  
# Set empty rectors for point-wise WAIC calculations
lppdv <- vector(length = n)
pwaicv <- vector(length = n)
  
# Loop through all counties
for (j in 1:n){
  Yi <- Y[j]
  mu <- rbeta %*% X[j, ]
  liki <- dnorm(Yi, mu, sqrt(rsig))
  lppdv[j] <- 1/B*sum(liki)
  pwaicv[j]	<- (1/(B-1))*sum((log(liki) - mean(log(liki)))^2)
}
  
# Calculate WAIC
lppd	<- sum(log(lppdv))
pwaic	<- sum(pwaicv)
WAIC_bio	<- -2*lppd +2*pwaic

# Display results
rbMat <- apply(rbeta, 2, quantile, probs = c(0.5, 0.025, 0.975))
postp <- apply(rbeta > 0, 2, mean)
outMat <- cbind(t(rbMat), postp)
rownames(outMat) <- colnames(X)
colnames(outMat) <- c(rownames(rbMat), 'P(b > 0)')
data_frame_bio <- as.data.frame(outMat)
data_frame_bio$Variables <- c('Intercept' , 'Influenza Percent 65+',
                              'Influenza Percent Female', 'Percent Vaccinated',
                              'Percent Obese', 'Percent Mental Distress',
                              'Avg. Humidity', 'Avg. Temperature')
# Save
write.csv(data_frame_bio, file = "tmp/data_frame_bio.csv")

# Human Error & Administrative Burden Explanation
# Reset B, n, Y, X, K
# Set the number of samples
B     <- 10000 
# Define the number of observations
n     <- nrow(regress_dat)
# Define the outcome vector
Y     <- regress_dat$percent_asymp_flu_scale

# Specify the design matrix
X     <- model.matrix(percent_asymp_flu_scale ~ 
                        # healthcare access variables
                        rate_doctors + rate_icu_beds + 
                        # SES variables
                        median_hh_income_scale + percent_white_flu_scale +
                        # healthcare burden variables
                        percent_flu + percent_rsv + percent_covid +
                        # urbanicity variable
                        percent_rural_scale, 
                      data = regress_dat)

# Define the number of independent variables (new design matrix)
K <- ncol(X)

# Calculate regression inputs
bhat <- c(solve(t(X)%*%X)%*%(t(X)%*%Y))
SSY	<- t(Y - X%*%bhat)%*%(Y - X%*%bhat)
XtXi <- solve(t(X)%*%X)

# Create a matric to store betas
rbeta <- matrix(0, nrow = B, ncol = K)

# Sample sigma values
rsig <- rinvgamma(B, (n-K)/2, (1/2)*SSY)

# Calculate a Var/Cov matrix from each sigma and sample betas
for (i in 1:B){
  CovX <- rsig[i]*XtXi
  rbeta[i,]	<- c(rmvn(1, mu = bhat, sigma = CovX))
}

# Set empty rectors for point-wise WAIC calculations
lppdv <- vector(length = n)
pwaicv <- vector(length = n)

# Loop through all counties
for (j in 1:n){
  Yi <- Y[j]
  mu <- rbeta %*% X[j, ]
  liki <- dnorm(Yi, mu, sqrt(rsig))
  lppdv[j] <- 1/B*sum(liki)
  pwaicv[j]	<- (1/(B-1))*sum((log(liki) - mean(log(liki)))^2)
}

# Calculate WAIC
lppd	<- sum(log(lppdv))
pwaic	<- sum(pwaicv)
WAIC_adm	<- -2*lppd +2*pwaic

# Display results
rbMat <- apply(rbeta, 2, quantile, probs = c(0.5, 0.025, 0.975))
postp <- apply(rbeta > 0, 2, mean)
outMat <- cbind(t(rbMat), postp)
rownames(outMat) <- colnames(X)
colnames(outMat) <- c(rownames(rbMat), 'P(b > 0)')
data_frame_adm <- as.data.frame(outMat)
data_frame_adm$Variables <- c('Intercept' , 'Rate Doctors',
                              'Rate ICU Beds', 'Median HH Income',
                              'Influenza Percent White', 'Rate Influenza',
                              'Rate COVID-19', 'Rate RSV', 'Percent Rural')
# Save
write.csv(data_frame_adm, file = "tmp/data_frame_adm.csv")

# Which WAIC is higher?
WAIC_adm < WAIC_bio

# Combine and order data
bar_data <- data.frame(Model=c('Immune/Biology', 'Burden/SES'),
                       WAIC=as.numeric(c(WAIC_bio, WAIC_adm)))

# Save
write.csv(data_frame_adm, file = "tmp/bar_plot.csv")

#########################################
# Variable selection: Forward Selection #
#########################################

# Reset B, n, Y, X, K
# Set the number of samples
B     <- 10000 
# Define the number of observations
n     <- nrow(regress_dat)
# Define the outcome vector
Y     <- regress_dat$percent_asymp_flu_scale

# Specify the design matrix for the full model
X     <- model.matrix(percent_asymp_flu_scale ~ 
                        # immunity variables
                        percent_65_plus_flu_scale + percent_female_flu_scale + 
                        percent_vaccinated_scale + percent_obese_scale + 
                        percent_mental_distress_scale + 
                        # biological variables
                        avg_humidity + avg_temp + 
                        # healthcare access variables
                        rate_doctors + rate_icu_beds + 
                        # SES variables
                        median_hh_income_scale + percent_white_flu_scale +
                        # healthcare burden variables
                        percent_flu + percent_rsv + percent_covid +
                        # urbanicity variable
                        percent_rural_scale, 
                      data = regress_dat)

# Define the number of independent variables (new design matrix)
K <- ncol(X)

# Calculate regression inputs
bhat <- c(solve(t(X)%*%X)%*%(t(X)%*%Y))
SSY	<- t(Y - X%*%bhat)%*%(Y - X%*%bhat)
XtXi <- solve(t(X)%*%X)

# Create a matric to store betas
rbeta <- matrix(0, nrow = B, ncol = K)

# Sample sigma values
rsig <- rinvgamma(B, (n-K)/2, (1/2)*SSY)

# Calculate a Var/Cov matrix from each sigma and sample betas
for (i in 1:B){
  CovX <- rsig[i]*XtXi
  rbeta[i,]	<- c(rmvn(1, mu = bhat, sigma = CovX))
}

# Display results
rbMat <- apply(rbeta, 2, quantile, probs = c(0.5, 0.025, 0.975))
postp <- apply(rbeta > 0, 2, mean)
dist_50 <- abs(postp - 0.5) # Calculate distance from 0.50
outMat <- as.data.frame(cbind(t(rbMat), postp, dist_50))
outMat$Variables <- colnames(X)
# Determine variable addition order
outMat |> arrange(desc(dist_50)) |> dplyr::select(Variables)

# STEPWISE #

# Set the number of samples
B <- 10000 
# Define the number of observations
n <- nrow(regress_dat)
# Define the outcome vector
Y <- regress_dat$percent_asymp_flu_scale

# Specify the design matrix
X <- model.matrix(percent_asymp_flu_scale ~ 
                    percent_65_plus_flu_scale + percent_flu + percent_rural_scale + 
                    percent_white_flu_scale + rate_icu_beds + avg_temp + 
                    percent_mental_distress_scale + percent_vaccinated_scale + 
                    percent_rsv + median_hh_income_scale + percent_obese_scale + 
                    rate_doctors + percent_female_flu_scale + percent_covid + 
                    avg_humidity, 
                  data = regress_dat)

# Set an empty vector for WAIC values
all_waic <- vector(length = 20)

# Loop through all variable combinations, calculating WAIC
for (i in 2:16) {
  # Print a tracking variable
  print(i)
  
  # Re-specify the design matrix based on the variable combination
  X_new <- X[, c(1:i)]
  
  # Define the number of independent variables (new design matrix)
  K <- ncol(X_new)
  
  # Calculate regression inputs
  bhat <- c(solve(t(X_new)%*%X_new)%*%(t(X_new)%*%Y))
  SSY	<- t(Y - X_new%*%bhat)%*%(Y - X_new%*%bhat)
  XtXi <- solve(t(X_new)%*%X_new)
  
  # Create a matric to store betas
  rbeta <- matrix(0, nrow = B, ncol = K)
  
  # Sample sigma values
  rsig <- rinvgamma(B, (n-K)/2, (1/2)*SSY)
  
  # Calculate a Var/Cov matrix from each sigma and sample betas
  for (j in 1:B){
    CovX <- rsig[j]*XtXi
    rbeta[j,]	<- c(rmvn(1, mu = bhat, sigma = CovX))
  }
  
  # Set empty rectors for point-wise WAIC calculations
  lppdv_test <- vector(length = n)
  pwaicv_test <- vector(length = n)
  
  # Loop through all counties
  for(h in 1:n){
    Yi <- Y[h]
    mu <- rbeta %*% X_new[h, ]
    liki <- dnorm(Yi, mu, sqrt(rsig))
    lppdv_test[h] <- 1/B*sum(liki)
    pwaicv_test[h]	<- (1/(B-1))*sum((log(liki) - mean(log(liki)))^2)
  }
  
  # Calculate WAIC
  lppd_test	<- sum(log(lppdv_test))
  pwaic_test	<- sum(pwaicv_test)
  WAIC	<- -2*lppd_test +2*pwaic_test
  
  # Add WAIC value to vector
  all_waic[i] <- WAIC
  
}

# Convert WAICs to data frame
WAIC <- as.data.frame(all_waic)
WAIC <- WAIC |> filter(all_waic != 0)

# Calculate pecent change
WAIC %>%
  mutate(pct_change = (all_waic/lag(all_waic) - 1) * 100)

# Save
write.csv(WAIC, file = "tmp/waic_forward.csv")
  
# Re-specify the design matrix based on the variable combination
X_new <- X[, c(1:9)]

# Define the number of independent variables (new design matrix)
K <- ncol(X_new)

# Calculate regression inputs
bhat <- c(solve(t(X_new)%*%X_new)%*%(t(X_new)%*%Y))
SSY	<- t(Y - X_new%*%bhat)%*%(Y - X_new%*%bhat)
XtXi <- solve(t(X_new)%*%X_new)

# Create a matric to store betas
rbeta <- matrix(0, nrow = B, ncol = K)

# Sample sigma values
rsig <- rinvgamma(B, (n-K)/2, (1/2)*SSY)

# Calculate a Var/Cov matrix from each sigma and sample betas
for (j in 1:B){
  CovX <- rsig[j]*XtXi
  rbeta[j,]	<- c(rmvn(1, mu = bhat, sigma = CovX))
}

# Display results
rbMat <- apply(rbeta, 2, quantile, probs = c(0.5, 0.025, 0.975))
postp <- apply(rbeta > 0, 2, mean)
outMat <- cbind(t(rbMat), postp)
rownames(outMat) <- colnames(X_new)
colnames(outMat) <- c(rownames(rbMat), 'P(b > 0)')
data_frame_for <- as.data.frame(outMat)
data_frame_for$Variables <- c('Intercept' , 'Influenza Percent 65+',
                              'Rate Influenza', 'Percent Rural',
                              'Influenza Percent White', 'Rate ICU Beds',
                              'Avg. Temperature', 'Percent Mental Distress',
                              'Percent Vaccinated')
# Save
write.csv(data_frame_for, file = "tmp/data_frame_for.csv")

################################################################################
################################################################################
