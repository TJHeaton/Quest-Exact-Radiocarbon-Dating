# Show density of annual tree-ring data over the Holocene
# Calendar ages are plotted in astronomical years 
# Overlain are the various Miyake Events identified in Holocene so far
# Magenta - Confirmed ESPEs
# Dark Green - Proposed ESPEs

# Store oldpar so do not change after running in main environment 
oldpar <- par(no.readonly = TRUE)

if(! "rintcal" %in% installed.packages())
  install.packages("rintcal")
require(rintcal)

dat <- read.table(system.file("extdata/", "intcal20_data.txt", package = "rintcal"), header=TRUE)
trees <- c(1:8, 59:99)
dat <- dat[which(dat$set %in% trees),] # trees only

yearly <- which(dat$calsig == 1)
yearlydata <- dat[yearly,]

DarkOrange <- rgb(213,94,0, (0.22*255), maxColorValue = 255)
DarkGreen <- rgb(0,158,115, (0.6*255), maxColorValue = 255)
Yellow <- rgb(240,228,66, (0.6*255), maxColorValue = 255)
Magenta <- rgb(1,0,1,.6) 

scalefac <- 1.2
ticksize <- 0.7


par(mar = c(3, 0.2, 0.1, 0.1))
plot(yearlydata$cal, yearlydata$calsig, 
     axes = FALSE,
     xlim = 1950 - c(13910,0), 
     ylim = c(0,1),
     type = "n",
     xlab = "Astronomical Year", 
     ylab = "", 
     bty= "l", 
     yaxs="i",
     mar=c(2,0,.5,.5), 
     mgp=c(2, .7, 0))

axis(1)
xtick <- seq(-14000, 2000, by = 200)
axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = 2/3 * -0.035)

# Miyake events:
known <- -c((7176 - 1), (5258 - 1), (660 - 1), -774, -993) 
poss <-  -c((5410 - 1), -1052, -1279)
segments(x0=known, 
         y0 = 0,
         y1 = 1.1,
         lwd = 5, 
         lend = 1,
         col= Magenta)
segments(x0 = poss,
         y0 = 0,
         y1 = 1.1,
         lwd = 5, 
         lend = 1,
         col= DarkGreen)


rug(1950 - yearlydata$cal, ticksize = ticksize, col=DarkOrange, lwd = 0.15)

# Revert to old par
par(oldpar)

