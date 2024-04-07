# Plot detail of IntCal20 vs IntCal13 over Thera period

# Include underlying tree-ring data
# X-axis (calendar age) bars represent the annual rings in the samples
# Y-axis (14C) bars represent the measurement uncertainty

# Store oldpar so do not change after running in main environment 
oldpar <- par(no.readonly = TRUE)

# Change par to allow extra space on LHS
par(mar = c(5, 4.4, 4, 2) + 0.1)

MookConvention <- TRUE
Mult_SigmaInterval <- 1 # Do we want 1 sigma or two sigma
y_lab_line <- 3.2 # Ensure does not overlap with labelling

# Select calendar age period to plot
range_calBC <- c(1720, 1480)
range_calBP <- 1950.5 + range_calBC

# Read in IntCal13 and IntCal20
IntCal13 <- read.csv("CalibrationCurves/intcal13.csv", header = TRUE, sep = ",") 
IntCal20 <- read.csv("CalibrationCurves/intcal20.csv", header = TRUE)

# Read in data - after repeat measurements have been merged
Extract <- read.csv("Datasets/TreeRings_IntCal20_MergedRepeats.csv", 
                    header = TRUE)
Extract$delta <- 1000*(Extract$F*exp(Extract$calage/8267) - 1)
Extract$deltaUpper <- 1000*((Extract$F + Mult_SigmaInterval *Extract$sdF) * exp(Extract$calage/8267) - 1)
Extract$deltaLower <- 1000*((Extract$F - Mult_SigmaInterval *Extract$sdF) * exp(Extract$calage/8267) - 1)
Extract$c14ageUpper <- Extract$c14age + Mult_SigmaInterval * Extract$c14sig
Extract$c14ageLower <- Extract$c14age - Mult_SigmaInterval * Extract$c14sig

# Identify samples to plot 
Use_IntCal20_ID <- which((IntCal20$CAL.BP < max(range_calBP)) & 
                           (IntCal20$CAL.BP > min(range_calBP)))
Selected_IntCal20 <- IntCal20[Use_IntCal20_ID,]
Selected_IntCal20$BCAge <- -(1950.5 - Selected_IntCal20$CAL.BP) 

plot(Selected_IntCal20$BCAge, Selected_IntCal20$c14_age, type = "n", col = "blue", 
     xlab = "Calendar Age (cal BC)", ylab = "",
     ylim = c(3140, 3510), xlim = range_calBC, 
     xaxs = "i", yaxs= "i",
     las = 1)

if(!MookConvention) {
  mtext(text = expression(paste("Radiocarbon Age (", ""^14, "C yr", " BP)")), 
        side = 2, line = y_lab_line) 
} else {
  mtext(text = expression(paste("Radiocarbon Age (BP)")), 
        side = 2, line = y_lab_line) 
}

# Add smaller ticks
xtick<-seq(1200, 1850, by=10)
ytick<-seq(3000, 4000, by=10)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)


# Plot shaded curve polygons before points
lines(-(1950.5 - IntCal13$calage), IntCal13$c14age, col = rgb(0, 102, 17, maxColorValue = 255))
polygon(-(1950.5 - c(IntCal13$calage, rev(IntCal13$calage))),
        c(IntCal13$c14age + Mult_SigmaInterval * IntCal13$c14sig, 
          rev(IntCal13$c14age -  Mult_SigmaInterval * IntCal13$c14sig)), 
        col= rgb(0, 102, 17, alpha = 127, maxColorValue = 255), 
        border=NA)

# IntCal20 over the top
lines(-(1950.5 - IntCal20$CAL.BP), IntCal20$X14C.age, col = rgb(213, 94, 0, maxColorValue = 255))
polygon(-(1950.5 - c(Selected_IntCal20$CAL.BP, rev(Selected_IntCal20$CAL.BP))),
        c(Selected_IntCal20$X14C.age + Mult_SigmaInterval * Selected_IntCal20$Sigma, 
          rev(Selected_IntCal20$X14C.age - Mult_SigmaInterval * Selected_IntCal20$Sigma)),
        col = rgb(213,94,0, alpha = 127, maxColorValue = 255), border=NA)

# Plot all obsverations points
points(-(1950.5 - Extract$calage), Extract$c14age, 
       pch = 19, cex = 0.2, 
       col = grey(level = 0, alpha = 0.9))
# Error-bars representing uncertainty in c14 axis
segments(x0 = -(1950.5 - Extract$calage), 
         y0 = Extract$c14ageUpper, 
         y1 = Extract$c14ageLower, 
         lwd = 0.3, col = grey(level = 0, alpha = 0.5))
# Bars representing the calendar years represented by the (multi-annual) blocked measurements
segments(x0 = -(1950.5 - floor(Extract$calage - (Extract$calspan-1)/2)), 
         x1 = -(1950.5 - floor(Extract$calage + (Extract$calspan-1)/2)),
         y0 = Extract$c14age, 
         lwd = 0.3, col = grey(level = 0, alpha = 0.5))

# (Redraw blocked measurements) to show the multi-annual IntCal20 tree-ring measurements in dark green
Blocked_ID <- which(Extract$calspan > 1)
Blocked_Extract <- Extract[Blocked_ID,]
points(-(1950.5 - Blocked_Extract$calage), 
       Blocked_Extract$c14age, 
       pch = 19, cex = 0.2, 
       col = rgb(0, 102, 17, maxColorValue = 255))
# Error-bars representing uncertainty in c14 axis
segments(x0 = -(1950.5 - Blocked_Extract$calage), 
         y0 = Blocked_Extract$c14ageUpper, 
         y1 = Blocked_Extract$c14ageLower, 
         lwd = 0.3, 
         col = rgb(0, 102, 17, maxColorValue = 255))
# Bars representing the calendar years represented by the (multi-annual) blocked measurements
segments(x0 = -(1950.5 - floor(Blocked_Extract$calage - (Blocked_Extract$calspan-1)/2)), 
         x1 = -(1950.5 - floor(Blocked_Extract$calage + (Blocked_Extract$calspan-1)/2)),
         y0 = Blocked_Extract$c14age, 
         lwd = 0.3, 
         col = rgb(0, 102, 17, maxColorValue = 255))


## Finally overlay some individual realisations
# Read in realisations
individual_realisations <- read.csv("Realisations/IntCal20_Realisations_c14age.csv")
# Select the colour palette
library(RColorBrewer)
n_samp <- ncol(individual_realisations) - 1 # First column is calendar age grid 
col_list <- brewer.pal(n_samp, "Set3") # Choose nsamp colours

for(i in 2:ncol(individual_realisations)) {
  lines(-(1950.5 - individual_realisations[,1]), 
        individual_realisations[,i], 
        col = col_list[i-1])
}

# Add legend 
legend("topright", legend = c("IntCal20", "IntCal13"), 
       fill = c(rgb(213, 94, 0, alpha = 127, maxColorValue = 255),
                rgb(0, 102, 17, alpha = 127, maxColorValue = 255)))



# Revert to old par
par(oldpar)







