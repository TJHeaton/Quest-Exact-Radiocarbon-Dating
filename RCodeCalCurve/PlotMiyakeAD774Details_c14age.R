# Plot detail of IntCal20 vs MiyakeCal over AD 744 period
# Plots in radiocarbon age (c14age) space

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

## Read in the data
Miyake <- 1950.5 - c(774.5)
# Read in data - after repeat measurements have been merged
Extract <- read.csv("Datasets/TreeRings_IntCal20_MergedRepeats.csv", 
                    header = TRUE)
Extract$delta <- 1000*(Extract$F*exp(Extract$calage/8267) - 1)
Extract$deltaUpper <- 1000*((Extract$F + Mult_SigmaInterval *Extract$sdF) * exp(Extract$calage/8267) - 1)
Extract$deltaLower <- 1000*((Extract$F - Mult_SigmaInterval *Extract$sdF) * exp(Extract$calage/8267) - 1)
Extract$c14ageUpper <- Extract$c14age + Mult_SigmaInterval * Extract$c14sig
Extract$c14ageLower <- Extract$c14age - Mult_SigmaInterval * Extract$c14sig

# Restrict plot to data 450 years either side of the 774AD Miyake Event 
Extract <- Extract[(Extract$calage < (Miyake + 450)) &  (Extract$calage > (Miyake - 450)),]
t_grid <- seq(min(Extract$calage),max(Extract$calage), by = 1)

# Read in IntCal20 curve and reduce to only section needed to plot
IntCal20 <- read.csv("CalibrationCurves/intcal20.csv", header = TRUE)
MiyakeID <- which((IntCal20$CAL.BP < Miyake + 70) & (IntCal20$CAL.BP > Miyake - 70))
Miyake20 <- IntCal20[MiyakeID,]
Miyake20$ADAge <- 1950.5 - Miyake20$CAL.BP 

# Create plot of curves
plot(Miyake20$ADAge, Miyake20$Delta.14C, type = "n", col = "blue", 
     xlab = "Calendar Age (cal AD)", ylab = "",
     ylim = c(1070, 1370), 
     xlim = 774.5 + c(-50,50),
     yaxs = "i",
     las = 1)

if(!MookConvention) {
  mtext(text = expression(paste("Radiocarbon Age (", ""^14, "C yr", " BP)")), 
        side = 2, line = y_lab_line) 
} else {
  mtext(text = expression(paste("Radiocarbon Age (BP)")), 
        side = 2, line = y_lab_line) 
}

# Add smaller ticks
xtick<-seq(700, 900, by=2)
ytick<-seq(1000, 1400, by=10)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)

# Plot observations
points(1950.5 - Extract$calage, Extract$c14age, pch = 19, cex = 0.3, col = grey(level = 0, alpha = 0.7))
# Error-bars representing uncertainty in Delta14C axis
segments(x0 = 1950.5 - Extract$calage, 
         y0 = Extract$c14ageUpper, 
         y1 = Extract$c14ageLower, 
         lwd = 0.3, col = grey(level = 0, alpha = 0.5))
# Bars representing the calendar years represented by the (multi-annual) blocked measurements
segments(x0 = 1950.5 - floor(Extract$calage - (Extract$calspan-1)/2), 
         x1 = 1950.5 - floor(Extract$calage + (Extract$calspan-1)/2),
         y0 = Extract$c14age, 
         lwd = 0.3, col = grey(level = 0, alpha = 0.5))

# Overlay shaded IntCal20 curve
polygon(1950.5 - c(IntCal20$CAL.BP, rev(IntCal20$CAL.BP)),
        c(IntCal20$X14C.age + Mult_SigmaInterval * IntCal20$Sigma, 
          rev(IntCal20$X14C.age - Mult_SigmaInterval * IntCal20$Sigma)),
        col = rgb(213,94,0, alpha = 127, maxColorValue = 255), border=NA)




# Overlay MiyakeCal curve
MiyakeCal <- read.csv("CalibrationCurves/MiyakeCalAD774Curve.csv")

MiyakeCal$c14age_upper <- MiyakeCal$c14age + Mult_SigmaInterval * MiyakeCal$c14sig
MiyakeCal$c14age_lower <- MiyakeCal$c14age - Mult_SigmaInterval * MiyakeCal$c14sig

polygon(1950.5 - c(MiyakeCal$calage, rev(MiyakeCal$calage)),
        c(MiyakeCal$c14age_upper, 
          rev(MiyakeCal$c14age_lower)),
        col=rgb(0,0,1,.5), border=NA)


## Finally overlay some individual MiyakeCal realisations
# Read in realisations
individual_realisations <- read.csv("Realisations/MiyakeCalAD774_Realisations_c14age.csv")
# Select the colour palette
library(RColorBrewer)
n_samp <- ncol(individual_realisations) - 1 # First column is calendar age grid 
col_list <- brewer.pal(n_samp, "Set3") # Choose nsamp colours

for(i in 2:ncol(individual_realisations)) {
  lines(1950.5 - individual_realisations[,1], 
        individual_realisations[,i], 
        col = col_list[i-1])
}

# Overlay mean delta14c curve
lines(1950.5 - MiyakeCal$calage, MiyakeCal$c14age, col = "blue")

legend("topright", legend = c("MiyakeCal", "IntCal20"), 
       fill = c(rgb(0,0,1,.5), 
                rgb(213,94,0, alpha = 127, maxColorValue = 255)))


# Revert to old par
par(oldpar)




