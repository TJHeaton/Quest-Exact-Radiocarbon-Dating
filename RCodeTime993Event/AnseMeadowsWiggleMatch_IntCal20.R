# Wiggle matching L'Anse aux Meadows Beams vs IntCal20 curve (which includes both Japenese and Buntgen reference data)

# Posterior for wiggle match is in red
# rgb(1,0,0,.5)
# NOTE: AGAIN THIS IS DARKER WITH alpha = 0.5 than Figures 1 and 2

AD <- TRUE
MookConvention <- TRUE
set.seed(17)


Mult_SigmaInterval <- 1 # Do we want 1 sigma or two sigma

# Wiggle matching L'Anse aux Meadows over a MiyakeCal curve using IntCal20 data
MiyakeCal <- read.csv("CalibrationCurves/intcal20.csv", header = TRUE)
names(MiyakeCal) <- c("calage", "c14age", "c14sig", "NA1", "NA2")
MiyakeCal$C14Upper <- MiyakeCal$c14age + Mult_SigmaInterval * MiyakeCal$c14sig
MiyakeCal$C14Lower <- MiyakeCal$c14age - Mult_SigmaInterval * MiyakeCal$c14sig
MiyakeCal$AD <- (1950 - MiyakeCal$calage) 

# Calibrate vs MiyakeCal
calcurve <- MiyakeCal 


beam <- "4A_68_E"

if(beam == "4A_59") { 
  wiggle_data <- read.csv("Anse_aux_Meadows/Beam4A_59_Meadows_WiggleMatch.csv")
} else if (beam == "4A_68_E") {
  wiggle_data <- read.csv("Anse_aux_Meadows/Beam4A_68_E_Meadows_WiggleMatch.csv")
} else if (beam == "4A_68_J") {
  wiggle_data <- read.csv("Anse_aux_Meadows/Beam4A_68_J_Meadows_WiggleMatch.csv")
} else {
  stop()
}

# Now create a function which calibrates multiple samples in a wiggle match
# Arguments:
# t - the age you want to test
# calgaps - the known claendar age gaps
# y - vector of c14ages
# er - observcational uncertainty
# calcurve - the calibration curve
Miyakecalibwiggle <- function(t, y, er, calgaps, calcurve) {
  # Find approx value of calcurve at
  wiggle_calages <- t + calgaps # ages of all samples
  wiggle_mean <- approx(calcurve[,1], calcurve[,2], xout = wiggle_calages, rule = 2)$y
  wiggle_sd <- approx(calcurve[,1], calcurve[,3], wiggle_calages, rule = 2)$y
  prod(dnorm(y, mean = wiggle_mean, sd = sqrt(er^2 + wiggle_sd^2)))                    
} 

fromto <- seq(900, 1000, by = 1)
ntheta <-length(fromto)
probs <- rep(NA, ntheta)
for(i in 1:ntheta) {
  probs[i] <- Miyakecalibwiggle(t = fromto[i], 
                                y = wiggle_data$c14_age, 
                                er = wiggle_data$c14_sig,
                                calgaps = wiggle_data$gaps,
                                calcurve = MiyakeCal)
}


# Normalise the posterior probabilities
probs <- probs/sum(probs)
maxcalprob <- max(probs)

# Check we get annual precision
probs[(which.max(probs) - 2):(which.max(probs) + 2)]
# which cal yr BP it is
fromto[which.max(probs)]

# Now make the plot
Miyake <- 1950 - 992.5 # c(774.5, 993.5)
ylimplot <- c(980, 1170)
xlimplot <- c(970, 1025)

plot(MiyakeCal$AD, MiyakeCal$c14age,
     xlab = "", ylab = "", 
     las = 1,  cex.lab = 1.3,
     xaxs="i", yaxs = "i", 
     type = "n",
     ylim=ylimplot , xlim = xlimplot)

mtext(text="Calendar Age (cal AD)", side=1, line=2, cex = cex_lab)
mtext(text = expression(paste("Radiocarbon Age (BP)")), 
      side = 2, line = y_lab_line, cex = cex_lab)


xtick<-seq(900, 1100, by=1)
ytick<-seq(900, 1400, by=10)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)

# Plot complete IntCal curve

lines(MiyakeCal$AD, MiyakeCal$c14age, 
      col = rgb(213,94,0, max = 255))
polygon(c(MiyakeCal$AD, rev(MiyakeCal$AD)), 
        c(MiyakeCal$C14Upper, rev(MiyakeCal$C14Lower)), 
        col = rgb(213,94,0, max = 255, alpha = 127), border=NA)

# Plot observations at their optimal locations
post_cal_ages_AD <- 1950- (fromto[which.max(probs)] + wiggle_data$gaps)
points(post_cal_ages_AD, wiggle_data$c14_age,  pch = 19, cex = 0.8)
segments(x0 = post_cal_ages_AD, 
         y0 = wiggle_data$c14_age + Mult_SigmaInterval *wiggle_data$c14_sig,
         y1 = wiggle_data$c14_age - Mult_SigmaInterval *wiggle_data$c14_sig)



# Create a histogram of the calibrated ages
pol <- cbind(fromto, probs)

# Manually adjust as plot function interpolates from adjacent years 
# Achieve by padding with 0 at every half year 
padprobs <- as.vector(rbind(0,v=probs))
padfromto <- as.vector(rbind(fromto - 0.5, fromto))
pol <- cbind(padfromto, padprobs)


if(AD) {
  polBC <- pol
  polBC[,1] <- (1950 - polBC[,1])
  pol <- polBC
}
pol[,1] <- pol[,1]
# Rescale to fit on the graph
pol[,2] <- pol[,2] * 0.8 * (ylimplot[2] - ylimplot[1]) / maxcalprob
pol[,2] <- pol[,2] + ylimplot[1]
polygon(pol, col = rgb(1,0,0,.5))

legend("topright", legend = c("IntCal20"), 
       lty = c(-1), pch = c(15), lwd = 2,
       col = c(rgb(213,94,0, max = 255, alpha = 127)), cex = 0.8, pt.cex = pt_cex)

# In this case we don't get annual precision (i.e., not greater than 95% posterior probability in a single year)
probs[(which.max(probs) - 2):(which.max(probs) + 2)] 
(1950 - fromto[(which.max(probs) - 2):(which.max(probs) + 2)])
# In this case 46.9%  in 1022 cal AD and 48.3% 1021 cal AD

text(x = 1021, y = max(pol[,2]) - 5, 
     labels = " 48% 1021 cal AD\n 47% 1022 cal AD", 
     cex = 0.8, pos = 2)















