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

## Variable selection: All combinations ##

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
all_waic <- vector(length = nrow(combo_list))

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
    liki <- dnorm(Yi, mu, rsig)
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
min(all_waic)

write.csv(as.data.frame(all_waic), file = "/Users/rcorgel/Downloads/WAIC.csv")

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
  liki <- dnorm(Yi, mu, rsig)
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

colors <- c('0-4' = '#c36272', '5-12' = '#C7A939', '13-17' = '#5a9374', 
            '18-49' = '#347dab', '50-64' = '#9086ba', '65+' = 'darkgrey'
)

bio_plot <- ggplot(data=data_frame_bio[data_frame_bio$Variables != 'Intercept', ], 
       aes(x=Variables, y= `50%`, ymin= `2.5%`, ymax= `97.5%`)) +
  geom_pointrange(aes(), color = '#5a9374', size = 1, linewidth = 1) + 
  ggtitle('Immune/Biology Model') + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Median (2.5%-97.5%)") +
  ylim(-1.6, 1.6) + 
  theme_minimal() +  # use a white background
  theme(legend.position = 'bottom') +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=28, hjust = 0.5))

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
  liki <- dnorm(Yi, mu, rsig)
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

colors <- c('0-4' = '#c36272', '5-12' = '#C7A939', '13-17' = '#5a9374', 
            '18-49' = '#347dab', '50-64' = '#9086ba', '65+' = 'darkgrey'
)

adm_plot <- ggplot(data=data_frame_adm[data_frame_adm$Variables != 'Intercept', ], 
       aes(x=Variables, y= `50%`, ymin= `2.5%`, ymax= `97.5%`)) +
  geom_pointrange(aes(), color = '#9086ba', size = 1, linewidth = 1) + 
  ggtitle('Burden/SES Model') + 
  ylim(-1.6, 1.6) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Median (2.5%-97.5%)") +
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=28, hjust = 0.5))

# Which WAIC is higher?
WAIC_adm < WAIC_bio

# Combine and order data
bar_data <- data.frame(Model=c('Immune/Biology', 'Burden/SES'),
                       WAIC=as.numeric(c(WAIC_bio, WAIC_adm)))

# Plot
bar_plot <- ggplot(data=bar_data, aes(x=Model, y=WAIC, fill = Model)) + 
  geom_bar(stat="identity", position="identity", color = 'black', 
           width = 0.5, alpha = 1) + 
  theme_classic() +
  scale_fill_manual('', values = c('#5a9374', '#9086ba'),
                    breaks = c('Immune/Biology', 'Burden/SES')) + 
  ylab('WAIC') + xlab("") + ggtitle('Model Comparison') +
  coord_cartesian(ylim = c(24000, 25000)) +
  theme(legend.text = element_text(size = 24),
                         legend.title = element_text(size = 20),
                         axis.text = element_text(size=20),
                         axis.title = element_text(size=24),
                         legend.position = 'none',
                         legend.box="vertical",
                         legend.margin=margin(),
                         strip.background = element_blank(),
                         legend.spacing.y = unit(0.25, 'cm'),
                         legend.key.size = unit(1, 'cm'),
                         strip.text = element_text(size = 16),
                         plot.title = element_text(size=28, hjust = 0.5))

figure_2 <- cowplot::plot_grid(adm_plot, bio_plot, ggplot + theme_void(), bar_plot,
                               nrow = 1, rel_widths = c(1, 1, 0.1, 1),
                               labels = c('', '', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_2.jpg', plot = figure_2, height = 5, width = 20)

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
                        percent_vaccinated_scale, 
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
    liki <- dnorm(Yi, mu, rsig)
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

WAIC <- as.data.frame(all_waic)
WAIC <- WAIC |> filter(all_waic != 0)

forward_line <- ggplot() + geom_line(aes(y = WAIC$all_waic, x = 1:15), color = '#c36272', 
                                     linewidth = 2) + 
  xlab('Step Number') + ggtitle('WAIC by Forward Selection Step') +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10, 12, 14)) +
  ylab('WAIC') + theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=28, hjust = 0.5))
  
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
                              
forward_plot <- ggplot(data=data_frame_for[data_frame_for$Variables != 'Intercept', ], 
                   aes(x=Variables, y= `50%`, ymin= `2.5%`, ymax= `97.5%`)) +
  geom_pointrange(aes(), color = '#c36272', size = 1, linewidth = 1) + 
  ggtitle('Forward Model') + 
  ylim(-1.6, 1.6) + 
  geom_hline(yintercept=0, lty=2) +  # add a dotted line at x=1 after flip
  coord_flip() +  # flip coordinates (puts labels on y axis)
  xlab("") + ylab("Median (2.5%-97.5%)") +
  theme_minimal() +  # use a white background
  theme(legend.position = 'bottom') +
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=28, hjust = 0.5))

figure_3 <- cowplot::plot_grid(forward_line, forward_plot,
                               nrow = 1, rel_widths = c(1, 1),
                               labels = c('', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_3.jpg', plot = figure_3, height = 5, width = 20)


library(sf)
# Load maps
usa_albers_state <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_state.rds'))   # convert to sf
usa_albers_county <- st_as_sf(readRDS('/Users/rcorgel/Library/CloudStorage/GoogleDrive-rcc92@georgetown.edu/.shortcut-targets-by-id/1Iyvoddzrygu8ZOPlcXa4H5Zh35AMBnZD/Ronan_Bansal_Lab/Projects/syndromic-surveillance-proj/tmp/usa_albers_county.rds')) # convert to sf

usa_albers_county <- left_join(usa_albers_county, regress_dat, by = c('GEOID' ='county_fips'))

map <- ggplot() +
  geom_sf(data = usa_albers_county, aes(fill = percent_asymp_flu_scale/100, group = GEOID), 
          color= 'black', linewidth = 0.15) +
  geom_sf(data = usa_albers_state, aes(group = STATEFP), 
          fill = '#FFFFFF00', color= 'black', linewidth = 0.4) +
  scale_fill_gradient('Percent\n', low = "white", high = "#347dab") + 
  ggtitle('Spatial Patterns') + 
  theme_void() + theme(legend.position = 'right',
                       legend.text = element_text(size = 20),
                       legend.title = element_text(size = 20),
                       axis.title = element_text(size=24),
                       legend.margin=margin(),
                       strip.background = element_blank(),
                       legend.spacing.y = unit(0.25, 'cm'),
                       legend.key.height = unit(1.5, 'cm'),
                       strip.text = element_text(size = 16),
                       plot.title = element_text(size=28, hjust = 0.5))

density <- ggplot(data = regress_dat) +
  geom_density(aes(percent_asymp_flu_scale /100), color = '#347dab',
               fill = '#347dab', alpha = 0.65) + 
  ggtitle('County Distribution') + 
  geom_vline(xintercept=mean(regress_dat$percent_asymp_flu_scale /100), 
             lty=2, linewidth = 1) + 
  ylab("Density") + xlab("Percent") + xlim(0, 1) + 
  theme_minimal() +  # use a white background
  theme(legend.text = element_text(size = 24),
        legend.title = element_text(size = 20),
        axis.text = element_text(size=20),
        axis.text.y = element_blank(),
        axis.title = element_text(size=24),
        legend.position = "bottom",
        legend.box="vertical",
        legend.margin=margin(),
        strip.background = element_blank(),
        legend.spacing.y = unit(0.25, 'cm'),
        legend.key.size = unit(1, 'cm'),
        strip.text = element_text(size = 16),
        plot.title = element_text(size=28, hjust = 0.5))

figure_1 <- cowplot::plot_grid(map, density,
                               nrow = 1, rel_widths = c(1, 1),
                               labels = c('', ''),
                               label_size = 26, hjust = 0) 

ggsave('./figs/figure_1.jpg', plot = figure_1, height = 6, width = 20)

################################################################################
################################################################################
