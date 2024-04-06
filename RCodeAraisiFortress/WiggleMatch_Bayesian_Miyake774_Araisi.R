# 28th April 2023
# Create plots for Nature perspective to show
# variable precision based on calibration curve

################################
# Plots required are:
# a) Calibration of single date of thera
# b) Calibration over Hallstat plateau (two 14C dates with same precision)
# c) Wiggle match over Hallstat plateau

# For this plot, I have used
# Curve rgb(0,0,1,.5)
# NOTE: IT IS DARKER THAN INTCAL20 FOR FIGS 1 and 2 (have chosen alpha = 0.5)

# Posterior for wiggle match is in red
# rgb(1,0,0,.5)
# NOTE: AGAIN THIS IS DARKER WITH alpha = 0.5 than Figures 1 and 2

# If you need another colour (for e.g. calibration of individual samples in animation) 
# then palette suggests 
# Yellow rgb(240,228,66, (0.2*255), maxColorValue = 255)


#library(fda)
#library("rstan")
#library(mvtnorm) # To sample from a mvn for the noise on t 
#library(Matrix)

AD <- TRUE
MookConvention <- TRUE
set.seed(17)

Mult_SigmaInterval <- 1 # Do we want 1 sigma or two sigma

# Wiggle matching Araisi Lake Fortress over a MiyakeCal curve
MiyakeCal <- read.csv("CalibrationCurves/MiyakeCalCurve.csv", header = TRUE)
MiyakeCal$C14Upper <- MiyakeCal$c14age + Mult_SigmaInterval * MiyakeCal$c14sig
MiyakeCal$C14Lower <- MiyakeCal$c14age - Mult_SigmaInterval * MiyakeCal$c14sig
MiyakeCal$AD <- (1950 - MiyakeCal$calage) 

# Calibrate vs MiyakeCal
calcurve <- MiyakeCal 

wiggle_data <- read.csv("AraisiWiggleMatch.csv")

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

fromto <- seq(1100, 1200, by = 1)
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
Miyake <- 1950 - 774.5 # c(774.5, 993.5)
ylimplot <- c(1088, 1370)
xlimplot <- 774.5 + c(-50,50)

par(mar = c(3, 4.4, 0.5, 0.5) + 0.1)
plot(MiyakeCal$AD, MiyakeCal$c14age,
     xlab = "", ylab = "", 
     las = 1,  cex.lab = 1.3,
     xaxs="i", yaxs = "i", 
     type = "n",
     ylim=ylimplot , xlim = xlimplot)

mtext(text="Calendar Age (cal AD)", side=1, line=2, cex = cex_lab)
mtext(text = expression(paste("Radiocarbon Age (BP)")), 
      side = 2, line = y_lab_line, cex = cex_lab)


xtick<-seq(700, 900, by=1)
ytick<-seq(1000, 1400, by=10)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)
mtext(LETTERS[4], side = 3, adj = labadj, cex = 4/3, font = 2,
      line = 4/3 * label_adj, outer = FALSE)


lines(MiyakeCal$AD, MiyakeCal$c14age, col = "blue")
polygon(c(MiyakeCal$AD, rev(MiyakeCal$AD)), 
        c(MiyakeCal$C14Upper, rev(MiyakeCal$C14Lower)), 
        col = rgb(0,0,1, 0.5), border=NA)
     

# Plot observations at their optimal locations
post_cal_ages_AD <- 1950 - (fromto[which.max(probs)] + wiggle_data$gaps)
points(post_cal_ages_AD, wiggle_data$c14_age,  pch = 19, cex = 0.8)
segments(x0 = post_cal_ages_AD, 
         y0 = wiggle_data$c14_age + Mult_SigmaInterval *wiggle_data$c14_sig,
         y1 = wiggle_data$c14_age - Mult_SigmaInterval *wiggle_data$c14_sig)



# Create a histogram of the calibrated ages
pol <- cbind(fromto, probs)

# Manually adjust as it interpolates from adjacent years 
# This falsely suggests some probability in those years too
max_index <- which.max(probs)
pol <- rbind(c(min(fromto), 0),
             pol[1:(max_index - 1),],
             c(pol[max_index,1] - 0.5, 0),
             pol[max_index,], 
             c(pol[max_index,1] + 0.5, 0),
             pol[(max_index+1):(dim(pol)[1]),], 
             c(max(fromto), 0))


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

legend("topright", legend = c("MiyakeCal"), 
       lty = c(-1), pch = c(15), lwd = 2,
       col = c(rgb(0,0,1,.5)), cex = 4/3, pt.cex = 2)

# mtext(paste0("(", LETTERS[3], ")"), side = 3, adj = 0.05, cex = 1.3, 
#       line = -1.3)

text(x = 782, y = max(pol[,2]), 
     labels = "99.9% 782 cal AD", 
     cex = 4/3, pos = 4)

# Adjust plot to create subplot 
u <- c(738, 770, 1120, 1245)
v <- c(
  grconvertX(u[1:2], "user", "ndc"),
  grconvertY(u[3:4], "user", "ndc")
)
par(fig=v, new=TRUE, mar=c(0,0,0,0) )

plot(MiyakeCal$AD, MiyakeCal$c14age,
     type = "l",
     col = "blue",
     las = 1,
     xlim = c(768, 780), 
     ylim = c(1150, 1320))
polygon(c(MiyakeCal$AD, rev(MiyakeCal$AD)), 
        c(MiyakeCal$C14Upper, rev(MiyakeCal$C14Lower)), 
        col = rgb(0,0,1, 0.5), border=NA)
xtick<-seq(700, 900, by=1)
ytick<-seq(1000, 1400, by=10)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)

# Put on shading in background
tempcal <- c(774,775) 
tempx <- c(tempcal, rev(tempcal))
tempy <- rep(par("usr")[3:4], c(2,2))
polygon(tempx, tempy, 
        border = NA,
        col = rgb(86, 180, 233, max = 255, alpha = 75))












