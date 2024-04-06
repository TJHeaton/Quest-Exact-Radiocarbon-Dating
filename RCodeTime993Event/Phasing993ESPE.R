# Case Study II - Regional differences in timing of 993 AD event
# and calibration of Anse aux Meadows
# Shows NH1, Japan, and IntCal20

# Store oldpar so do not change after running in main environment 
oldpar <- par(no.readonly = TRUE)

# Plotting parameters
labadj <- 0.01
cex_lab <- 0.8
y_lab_line <- 3.2
scalefac <- 1.3
label_adj <- -1.5
pt_cex <- 0.8

# Read in Buntgen and Miyake Data 
Buntgen_NH0 <- read.csv("ObservationalData/Buntgen2018_NH0_Data.csv", header = TRUE)
Buntgen_NH1 <- read.csv("ObservationalData/Buntgen2018_NH1_Data.csv", header = TRUE)
Buntgen_NH2 <- read.csv("ObservationalData/Buntgen2018_NH2_Data.csv", header = TRUE)
Miyake_Japan <- read.csv("ObservationalData/Miyake_993_Data.csv", header = TRUE)

# Calculate Delta14C for Buntgen and plot
Buntgen_NH0$Delta14C <- 1000 * (exp(Buntgen_NH0$calage/8267) * exp(-Buntgen_NH0$c14age/8033) - 1)
Buntgen_NH1$Delta14C <- 1000 * (exp(Buntgen_NH1$calage/8267) * exp(-Buntgen_NH1$c14age/8033) - 1)
Buntgen_NH2$Delta14C <- 1000 * (exp(Buntgen_NH2$calage/8267) * exp(-Buntgen_NH2$c14age/8033) - 1)

Buntgen_NH0$Delta14C_Sig <- Buntgen_NH0$c14sig/8.033
Buntgen_NH1$Delta14C_Sig <- Buntgen_NH1$c14sig/8.033
Buntgen_NH2$Delta14C_Sig <- Buntgen_NH2$c14sig/8.033

Buntgen_NH0$Zone <- "NH0 - Buntgen et al. (2018)"
Buntgen_NH1$Zone <- "NH1 - Buntgen et al. (2018)"
Buntgen_NH2$Zone <- "NH2 - Buntgen et al. (2018)"
Miyake_Japan$Zone <- "Japan - Miyake et al. (2013,2014)"

# Function to combine measurements with same calendar age in each Zone
Combine_Measurements <- function(x) {
  # Create a dataframe to store the averaged values in
  TempData <- setNames(data.frame(matrix(ncol = 5, nrow = 0)), c("AD_age", "calage", "Delta14C", "Delta14C_Sig", "Zone"))
  calages <- unique(x$calage)
  zone <- unique(x$Zone)

  for(i in 1:length(calages)) {
    TempAge <- calages[i]
    AD_Age <- 1950 - TempAge
    AgeRep <- which(x$calage == TempAge)
    TempDelta14C <- x$Delta14C[AgeRep]
    TempDelta14C_Prec <- 1/x$Delta14C_Sig[AgeRep]^2
  
    # Choose the minimum unbiased estimator
    # Xbar = [(prec1 * X_1) + (prec2 * X_2) + ... + (precn * X_n)] / [ prec1 + prec2 + ... + precn]  
    Comb_Delta14C <- sum(TempDelta14C_Prec * TempDelta14C) / sum(TempDelta14C_Prec) # Average the Delta14C values
    Comb_Delta14C_Sig <- sqrt(1/sum(TempDelta14C_Prec)) # Add in quadrature
  
    # Now create new row storing this information (appending on to current output)
    CombData <- data.frame("AD_age" = AD_Age, 
                           "calage" = TempAge, 
                           "Delta14C" = Comb_Delta14C,
                           "Delta14C_Sig" = Comb_Delta14C_Sig, 
                           "Zone" = zone)
    TempData <- rbind(TempData, CombData)
  }
  return(TempData)
}

Combined_NH0 <- Combine_Measurements(Buntgen_NH0)
Combined_NH1 <- Combine_Measurements(Buntgen_NH1)
Combined_NH2 <- Combine_Measurements(Buntgen_NH2)
Combined_Miyake <- Combine_Measurements(Miyake_Japan)

# All sets combined 
Combined_993_Data <- rbind(Combined_NH0,
                           Combined_NH1,
                           Combined_NH2,
                           Combined_Miyake)

Combined_993_Data$Zone <- as.factor(Combined_993_Data$Zone)
Combined_993_Data$DeltaUpper <- Combined_993_Data$Delta14C + Combined_993_Data$Delta14C_Sig
Combined_993_Data$DeltaLower <- Combined_993_Data$Delta14C - Combined_993_Data$Delta14C_Sig

library(colorspace)
n_zones <- nlevels(Combined_993_Data$Zone)
col_list <- rainbow_hcl(n_zones)

NH0_col <- col_list[which(levels(Combined_993_Data$Zone) == "NH0 - Buntgen et al. (2018)")]
NH1_col <- col_list[which(levels(Combined_993_Data$Zone) == "NH1 - Buntgen et al. (2018)")]
Japan_col <- col_list[which(levels(Combined_993_Data$Zone) == "Japan - Miyake et al. (2013,2014)")]


layout(matrix(c(1,2, 1,3, 1, 4), 3, 2, byrow = TRUE))

par(mar = c(3, 4.4, 0.5, 0.5) + 0.1)

# Plot 1 - The different phasing of the 993 ESPE in different regions
plot(Combined_993_Data$AD_age, Combined_993_Data$Delta14C, type = "n", col = "blue", 
     xlab = "", ylab = "",
     ylim = c(-25,-5), las = 1)

# Put on shading in background
tempcal <- c(992,993) 
tempx <- c(tempcal, rev(tempcal))
tempy <- rep(par("usr")[3:4], c(2,2))
polygon(tempx, tempy, 
        border = NA,
        col = rgb(86, 180, 233, max = 255, alpha = 75))

tempcal <- c(993,994) 
tempx <- c(tempcal, rev(tempcal))
tempy <- rep(par("usr")[3:4], c(2,2))
polygon(tempx, tempy, 
        border = NA,
        col = rgb(240, 228, 66, max = 255, alpha = 75))

# Add axes labels
mtext(text="Calendar Age (cal AD)", side=1, line=2, cex = cex_lab)
mtext(text = expression(paste(Delta^14, "C", " (\u2030)")), 
      side = 2, line = y_lab_line - 0.2, cex = cex_lab)

xtick<-seq(980, 1010, by=1)
ytick<-seq(-30, 0, by=1)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.015)

# Now plot the data
points(Combined_993_Data$AD_age, Combined_993_Data$Delta14C, 
       pch = 19, cex = pt_cex, col = col_list[Combined_993_Data$Zone])
       
segments(x0 = Combined_993_Data$AD_age, 
         y0 = Combined_993_Data$DeltaUpper,
         y1 = Combined_993_Data$DeltaLower,
         lwd = 0.7, col = col_list[Combined_993_Data$Zone]) 
  
for(i in 1:n_zones) {
  zone_name <- levels(Combined_993_Data$Zone)[i]
  ZonalData <- subset(Combined_993_Data, Zone == zone_name)
  lines(ZonalData$AD_age, ZonalData$Delta14C, 
        col = col_list[i], lwd = 3)
}

legend("topright", 
       legend = levels(Combined_993_Data$Zone), 
       col = col_list, 
       pch = 19,
       cex = 0.8, 
       pt.cex = pt_cex)


# Now calibrate against various indiivdual datasets and plot 
source("AnseMeadowsWiggleMatch_Buntgen_NH1.R")
source("AnseMeadowsWiggleMatch_Japanese.R")
source("AnseMeadowsWiggleMatch_IntCal20.R")

# Revert to old par
par(oldpar)



