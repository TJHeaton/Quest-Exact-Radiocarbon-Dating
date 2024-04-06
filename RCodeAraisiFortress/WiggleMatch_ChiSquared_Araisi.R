# Perform a chi-squared wiggle-match of Lake Araisi timber fortress

# Read in Araisi data (from Meadows paper)
wiggle_data <- read.csv("AraisiData/AraisiWiggleMatch.csv")
# Reverse so that oldest dates are at the start
wiggle_data <- wiggle_data[order(wiggle_data$Year),]

# Read in Buntgen NH0 data
Buntgen_curve <- read.csv("CalibrationCurves/Buntgen774Event.csv", header = TRUE)
use_ID <- which(Buntgen_curve$AD_age >= 771 & Buntgen_curve$AD_age <= 779)
Buntgen_curve <- Buntgen_curve[use_ID,]
Buntgen_range <- range(Buntgen_curve$AD_age)
rm(use_ID)


n_araisi <- dim(wiggle_data)[1] 

# Now create 
wiggle_data$baseAD_age <- seq(from = min(Buntgen_range), by = 1, length = n_araisi)

chisquared_stat <- rep(NA, 6)
shifts <- -(0:5)
for(i in 1:6) {
  
  temp_shift <- shifts[i]
  shift_Araisi_age <- wiggle_data$baseAD_age + temp_shift
  # Find which calendar ages overlap with Buntgen range
  useID <- which(shift_Araisi_age >= Buntgen_range[1] & shift_Araisi_age <= Buntgen_range[2])
  araisi_overlap <- wiggle_data[useID, c("c14_age", "c14_sig")]
  if(dim(araisi_overlap)[1] != 9) stop("Error in overlap")
  
  # Now find the chi-squared stat 
  numerator <-  (araisi_overlap$c14_age - Buntgen_curve$c14age)^2 
  denominator <- araisi_overlap$c14_sig^2 + Buntgen_curve$c14sig^2
  chisquared_stat[i] <- sum(numerator/denominator)  
}

optimal_shift <- shifts[which.min(chisquared_stat)]
wiggle_data$optimal_AD_age <- wiggle_data$baseAD_age + optimal_shift

# Create a plot of the chi-squared statistic
plot(shifts, chisquared_stat,
     xlab = "", ylab = "", 
     las = 1,  cex.lab = 1.3,
     xaxs="i", yaxs = "i", 
     type = "n",
     main = expression(paste(chi^2, " wiggle match of Lake Āraiši Fortress")),
     ylim = c(0, 195),
     xlim = c(778.5, 784.5))

mtext(text="Year n - 53 (cal AD)", side=1, line=2, cex = cex_lab)
mtext(text = expression(paste(chi^2, "-statistic")), 
      side = 2, line = y_lab_line -1 , cex = cex_lab)

lines(784+shifts, chisquared_stat, col = "black", lwd = 3)
points(784+shifts, chisquared_stat, col = "black", pch = 19, cex = 1.5)
abline(h = qchisq(0.95, 8), col = "green", lty = 2, lwd = 3)

text(779.5, qchisq(0.95, 8), 
     labels = expression(paste(chi["0.95, 9"]^2, "= 15.5")), 
     pos = 3)




