# Function which calibrates multiple samples in an independent wiggle match
# Arguments:
# cal_age - the calendar age you want to test for sample with cal_gaps = 0
# radiocarbon_ages - vector of c14ages
# radiocarbon_sigmas - observational uncertainty
# cal_gaps - the known calendar age gaps between samples
# calcurve - the calibration curve
IndCalibWiggle <- function(cal_age, radiocarbon_ages, radiocarbon_sigmas, cal_gaps, calcurve) {
  wiggle_cal_ages <- cal_age + cal_gaps # calendar ages of all samples
  wiggle_mean <- approx(calcurve[,1], calcurve[,2], xout = wiggle_cal_ages, rule = 2)$y
  wiggle_sd <- approx(calcurve[,1], calcurve[,3], wiggle_cal_ages, rule = 2)$y
  posterior_likelihood <- prod(dnorm(radiocarbon_ages, 
                                     mean = wiggle_mean, 
                                     sd = sqrt(radiocarbon_sigmas^2 + wiggle_sd^2)))
  return(posterior_likelihood)
} 


# Function which calibrates multiple samples in a DEPENDENT wiggle match
# This calibration assumes that the data are offset from the cal curve mean
# Does not incorporate uncertainty on calibration curve
# Arguments:
# cal_age - the calendar age you want to test for sample with cal_gaps = 0
# radiocarbon_ages - vector of c14ages
# radiocarbon_sigmas - observational uncertainty (the independent component)
# cal_gaps - the known calendar age gaps between samples
# offset_sd - the uncertainty on the shared offset
# calcurve - the calibration curve
OffsetCalibWiggle <- function(cal_age, radiocarbon_ages, radiocarbon_sigmas, cal_gaps, offset_sd, calcurve) {
  n_obs <- length(radiocarbon_ages)
  wiggle_cal_ages <- cal_age + cal_gaps # calendar ages of all samples
  wiggle_mean <- approx(calcurve[,1], calcurve[,2], xout = wiggle_cal_ages, rule = 2)$y
  addcovar <- offset_sd^2 * matrix(1, nrow = n_obs, ncol = n_obs)
  posterior_log_likelihood <- dmvnorm(radiocarbon_ages, 
                    mean = wiggle_mean, 
                    sigma = diag(radiocarbon_sigmas^2) + addcovar, 
                    log = TRUE)
  posterior_likelihood <- exp(posterior_log_likelihood)
  return(posterior_likelihood)
} 
