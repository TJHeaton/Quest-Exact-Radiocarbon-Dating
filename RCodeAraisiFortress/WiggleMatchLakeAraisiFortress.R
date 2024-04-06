# Case Study 1 - Lake Araisi
# Plots showing three appproaches to calibrate of Lake Araisi Fortress

# Plot 1 (panel B in perspective) - IntCal20 calibration
# Plot 2 (panel C in perspective) - Chi-squared calibration
# Plot 3 (pabel D in perspective) - MiyakeCal calibration

labadj <- 0.01
cex_lab <- 0.8
y_lab_line <- 3.2
scalefac <- 1.3
label_adj <- -1

# Now panel B
source("PanelB_IntCal20_AraisiWiggleMatch_CBlind_1sigma.R")

source("ChiSquaredAraisi.R")

source("PanelD_Miyake774_AraisiWiggleMatch_CBlind_1sigma.R")
