

# OBTAIN KELT RATES AND ASSOCIATED 95% CONFIDENCE LIMITS VIA PARAMETRIC BOOTSTRAP
# BY: BEN STATON (CRITFC)
# FOR: LAURA JENKINS (U OF IDAHO)
# DATE: 9/25/2025

# set the random seed so the random number generator is the same every time this code is executed
set.seed(1234)

# enter data
dat = data.frame(
  year           = c(2016, 2017, 2018, 2019, 2020, 2021),
  total_tagged   = c(2798, 2117, 1403, 1094, 1448, 2275),
  total_detected = c( 138,   31,   78,   66,  275,  489),
  detection_mean = c(0.06, 0.02, 0.05, 0.07, 0.25, 0.42),
  detection_lwr  = c(0.04, 0.01, 0.03, 0.05, 0.20, 0.35),
  detection_upr  = c(0.09, 0.04, 0.07, 0.10, 0.30, 0.49)
)

# function to perform the logit transformation
# also known as "log-odds"; also given by qlogis() in R
logit = function(p) log(p/(1 - p))

# function to perform the inverse logit transformation
# also given by plogis() in R
ilogit = function(lp) exp(lp)/(1 + exp(lp))

# get SE(p) on logit-scale
# slightly different for upper and lower bounds, due to rounding when reporting the CIs of p
logit_p_se_from_upr = (logit(dat$detection_upr) - logit(dat$detection_mean))/1.96
logit_p_se_from_lwr = (logit(dat$detection_mean) - logit(dat$detection_lwr))/1.96

# take the mean of the two estimates by year to reconcile the difference
dat$logit_p_se = rowMeans(cbind(logit_p_se_from_lwr, logit_p_se_from_upr))

# simulate 1000 values of detection probability for each year based on the mean and SE
random_detections = lapply(1:nrow(dat), function(y) ilogit(rnorm(1000, logit(dat$detection_mean[y]), dat$logit_p_se[y])))

# perform the kelt rate calculation for each of the 1000 draws per year
random_kelt_rates = lapply(1:nrow(dat), function(y) {
  out = (dat$total_detected[y]/random_detections[[y]])/dat$total_tagged[y]
  
  # impose a maximum of 1: more individuals cannot kelt than were tagged
  # do it here so the frequency of 1's is properly accounted for in bootstrap summary
  out[out > 1] = 1
  
  # return output
  out
})

# summarize the bootstrap samples for each year
boot_summary = sapply(random_kelt_rates, function(x) c(boot_mean = mean(x), lwr95 = unname(quantile(x, 0.025)), upr95 = unname(quantile(x, 0.975))))

# obtain the point estimate based on no simulation: uses point estimate of detection rate rather than bootstrap distribution
data_estimates = (dat$total_detected/dat$detection_mean)/dat$total_tagged
data_estimates[data_estimates > 1] = 1

# combine with bootstrap summary
output = rbind(data_mean = data_estimates, boot_summary)
colnames(output) = dat$year

# convert to percentage
output = round(output, 4) * 100
output = apply(output, 2, function(x) paste0(x, "%"))

# format output for printing
output = as.data.frame(output)
rownames(output) = c("pt_ests_only", "boot_mean", "boot_lwr95", "boot_upr95")

# view the output
print(output)
