# Show effect of calibrating with and without shared offset

# Use colour-blind friendly palette from
# https://www.color-hex.com/color-palette/49436
# and plot 1 sigma intervals on the curve and the observations

# Illustrate the effect of an offset on 14C calibration over the Hallstatt Plateau
################################
# Approach required is:
# a) Create simulated data over Hallstatt Plateau by shifting IntCal20 mean
# b) Wiggle-match recognising that there is a shared offset (correct answer)
# c) Wiggle-match on basis uncertainties are independent (wrong answer) 


set.seed(17)
library(mvtnorm)
source("WiggleMatchFunctions.R")

# Store oldpar so do not change after running in main environment 
oldpar <- par(no.readonly = TRUE)

dens_plot_scale <- 0.2 # How far up you want the density plot to go
BC <- TRUE
MookConvention <- TRUE
Mult_SigmaInterval <- 1 # Do we want 1-sigma or 2-sigma intervals

# Number of samples in the wiggle match (likely more will illustrate better)
n_wiggle <- 10

# Parameters in sim study
offset <- -16 # Shared offset from mean IntCal20 curve
ind_sd <- 5 # Additional independent sd
cal_gap <- 5 # Gap between samples

# Plotting ranges
cal_age_lim_plot <- c(2340, 2750) 
radiocarbon_age_lim_plot <- c(2280.1, 2649.9)


# Calibrate against the IntCal20 curve
IntCal20 <- read.csv("intcal20.csv", header = TRUE)
IntCal20$C14Upper <- IntCal20$X14C.age + Mult_SigmaInterval * IntCal20$Sigma
IntCal20$C14Lower <- IntCal20$X14C.age - Mult_SigmaInterval * IntCal20$Sigma
IntCal20$BC <- -(1950.5 - IntCal20$CAL.BP) 
calcurve <- IntCal20 

# Create simulated data
truethetaBC <- 552.5
index <- which.min(abs(calcurve$BC - truethetaBC))
truethetacalBP <- calcurve$CAL.BP[index]
wiggle_calage <- seq(from = truethetacalBP, by = cal_gap, length.out = n_wiggle)
wiggle_c14mu <- approx(calcurve$CAL.BP, calcurve$X14C.age, xout = wiggle_calage)$y

# Create observed c14 ages that are offset (with extra bit of independent uncertainty)
wiggle_c14ages <- wiggle_c14mu + offset
wiggle_c14ages <- rnorm(n_wiggle, mean = wiggle_c14ages, sd = ind_sd)

# Create vector of known calendar age gaps
cal_gaps <- seq(from = 0, by = cal_gap, length.out = n_wiggle)

wiggle_data <- data.frame(gaps = cal_gaps, c14_age = wiggle_c14ages, c14_sig = ind_sd)



# Choose a suitable set of calendar ages to scan over 
min_cal_age <- 2000
max_cal_age <- 3000
cal_age_grid <- seq(min_cal_age, max_cal_age, by = 1)
n_theta <-length(cal_age_grid)

# Find posterior probabilities under independent calibration 
ind_posterior_probs <- rep(NA, n_theta)
for(i in 1:n_theta) {
  ind_posterior_probs[i] <- IndCalibWiggle(cal_age = cal_age_grid[i],
                             radiocarbon_ages = wiggle_data$c14_age, 
                             radiocarbon_sigmas = wiggle_data$c14_sig, 
                             cal_gaps = wiggle_data$gaps,
                             calcurve = calcurve)
}
# Normalise the posterior probabilities
ind_posterior_probs <- ind_posterior_probs/sum(ind_posterior_probs)
ind_max_posterior_prob <- max(ind_posterior_probs)


par(mar = c(5, 4.8, 2, 2) + 0.1)
plot(cal_age_grid, ind_posterior_probs, type="n", 
     xlab = "Calendar Age (cal BC)", ylab = "", 
     las = 1,  cex.lab = 1.3,
     xaxs ="i", yaxs = "i",
     ylim= radiocarbon_age_lim_plot , xlim = -(1950.5 - rev(cal_age_lim_plot)))
xtick <- seq(200,1000, by = 10)
axis(side = 1, at = xtick, labels = FALSE, lwd = 0.5, tck = -0.01)
ytick <- seq(2000, 3000, by = 10)
axis(side = 2, at = ytick, labels = FALSE, lwd = 0.5, tck = -0.01)
if(!MookConvention) {
  mtext(text = expression(paste("Radiocarbon Age (", ""^14, "C yr", " BP)")), side = 2, line = 3.2, cex = 1.3) 
} else {
  mtext(text = expression(paste("Radiocarbon Age (BP)")), side = 2, line = 3.2, cex = 1.3) 
}

# Plot calibration curve
lines(calcurve$BC,calcurve$X14C.age, 
      col = rgb(213,94,0, maxColorValue = 255))
polygon(c(calcurve$BC, rev(calcurve$BC)), 
        c(calcurve$C14Upper, rev(calcurve$C14Lower)), 
        col = rgb(213,94,0, alpha = 127, maxColorValue = 255),
        border=NA)

# True calendar age of most recent sample
abline(v = truethetaBC, lty = 2, lwd = 2, col = "red")


legend("topright", legend = c("IntCal20"), 
       lty = c(-1), pch = c(15), lwd = 2,
       col = c(rgb(213,94,0, alpha = 127, maxColorValue = 255)), cex = 1, pt.cex = 2)


# Now plot the posterior calendar age fit for the wiggle match

# Create a histogram of the calibrated ages
pol <- cbind(c(min(cal_age_grid), cal_age_grid, max(cal_age_grid)), c(0, ind_posterior_probs, 0))
if(BC) {
  polBC <- pol
  polBC[,1] <- -(1950.5 - polBC[,1])
  pol <- polBC
}
pol[,1] <- pol[,1] + min(wiggle_data$gaps)
# Rescale to fit on the graph
pol[,2] <- pol[,2] * dens_plot_scale * (radiocarbon_age_lim_plot[2] - radiocarbon_age_lim_plot[1]) / ind_max_posterior_prob
pol[,2] <- pol[,2] + radiocarbon_age_lim_plot[1]
polygon(pol, col = rgb(1,0,1,.5))


# Now repeat with offset wiggle match
offset_posterior_probs <- rep(NA, n_theta)
for(i in 1:n_theta) {
  offset_posterior_probs[i] <- OffsetCalibWiggle(cal_age = cal_age_grid[i], 
                                                 radiocarbon_ages = wiggle_data$c14_age,
                                                 radiocarbon_sigmas = wiggle_data$c14_sig, 
                                                 cal_gaps = wiggle_data$gaps,
                                                 offset_sd = 25, 
                                                 calcurve = calcurve)
}

# Normalise the posterior probabilities
offset_posterior_probs <- offset_posterior_probs/sum(offset_posterior_probs)

# Plot the observations on the curve
post_cal_ages <- cal_age_grid[which.max(offset_posterior_probs)] + cal_gaps
post_BC_ages <- -(1950.5 - post_cal_ages)
points(post_BC_ages, wiggle_c14ages,  pch = 19, cex  = 0.5)
segments(x0 = post_BC_ages, 
         y0 = wiggle_data$c14_age + Mult_SigmaInterval * wiggle_data$c14_sig, 
         y1 = wiggle_data$c14_age - Mult_SigmaInterval * wiggle_data$c14_sig)

## Now plot the posterior calendar age fit for the offset wiggle match

# Create a histogram of the calibrated ages
pol <- cbind(c(min(cal_age_grid), cal_age_grid, max(cal_age_grid)), c(0, offset_posterior_probs, 0))
if(BC) {
  polBC <- pol
  polBC[,1] <- -(1950.5 - polBC[,1])
  pol <- polBC
}
pol[,1] <- pol[,1] + min(wiggle_data$gaps)
# Rescale to fit on the graph
pol[,2] <- pol[,2] * dens_plot_scale * (radiocarbon_age_lim_plot[2] - radiocarbon_age_lim_plot[1]) / ind_max_posterior_prob
pol[,2] <- pol[,2] + radiocarbon_age_lim_plot[1]
polygon(pol, col = rgb(0,1,0,.5))

# Revert to old par
par(oldpar)






