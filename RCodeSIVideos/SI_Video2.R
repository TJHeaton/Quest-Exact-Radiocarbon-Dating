# load an R package which has the IntCal20 calibration curve
if(!'rintcal'%in% installed.packages())
  install.packages('rintcal')
require(rintcal) # but using Miyake cc

# can be found under RCodeCalCurve/CalibrationCurves/MiyakeCalAD774Curve.csv
cc <- read.table("../RCodeCalCurve/CalibrationCurves/MiyakeCalAD774Curve.csv", header=TRUE, sep=",")[,1:3]

cc[,1] <- 1950.5 - cc[,1]
cc.pol <- cbind(c(cc[,1], rev(cc[,1])), c(cc[,2]-cc[,3], rev(cc[,2]+cc[,3]))) 

# actual C-14 dates from Meadows et al. 2023 (Radiocarbon)
dat <- read.table("../RCodeAraisiFortress/AraisiData/AraisiWiggleMatch.csv")
# dat <- read.table("Araisi.dat")
dat <- dat[nrow(dat):1,] # reverse
n <- nrow(dat)
y <- dat[,1]
er <- dat[,2]
rings <- dat[,3] - min(dat[,3]) # most recent dated ring at 53 years; set to 0

# now calibrate
xlimplot <- c(700, 850)
xseq <- seq(xlimplot[1], xlimplot[2], length=1e3) # manual
ccy <- approx(cc[,1], cc[,2], xseq)$y # using Miyake cc
ccer <- approx(cc[,1], cc[,3], xseq)$y
cal <- array(0, dim=c(n, length(xseq)))
for(i in 1:n)
  cal[i,] <- dnorm(y[i], ccy, sqrt(er[i]^2 + ccer^2))

# find the probabilities belonging to a calendar age, based on the dates and their known spacing
find.wmd <- function(yr) {
  yrseq <- yr - rings # will change as yr changes
  wmd_mu <- approx(cc[,1], cc[,2], yrseq)$y
  wmd_sd <- approx(cc[,1], cc[,3], yrseq)$y
  return(prod(dnorm(y, wmd_mu, sqrt(wmd_sd^2+er^2))))
}
wmd.probs <- c()
for(i in 1:length(xseq))
  wmd.probs[i] <- find.wmd(xseq[i])
wmd.probs <- wmd.probs/max(wmd.probs) # set peak to 1

ex1 <- 40 # for the individual calibrated dates
ex2 <- 15 # for the 'wiggle match'
ex3 <- 300 # radiocarbon ages

xpol <- list()
for(i in 1:n)
  xpol[[i]] <- cbind(c(min(xseq), xseq, max(xseq)), i+ex1*c(0, cal[i,], 0))
yseq <- seq(min(y)-(5*max(er)), max(y)+(5*max(er)), length=1000)

plot.wmd <- function(yr, BCAD=TRUE, c14.lim=c()) {
  op <- par(bty="o", mar=c(3,5,1,1), mgp=c(2,.7,0), new=FALSE, xaxs="i", las=1)
  plot(0, type="n", xlim=range(xseq), xlab="Calendar Age (AD)", ylim=c(0, 35), ylab="", yaxt="n", yaxs="i")
  xtick <- seq(500, 900, by=2)
  axis(side=1, at=xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
  legend("topright", legend="MiyakeCal", fill=rgb(0,0,1,.5), border=NA)

  for(i in 1:n) {
    # plot the calibrated distributions
    polygon(xpol[[i]], col=rgb(0,.5,0,.25), border=rgb(0,.5,0,.25))
    abline(h=i, col=rgb(0,.5,0,.25))

    # draw lines and points to show where we are on the calibrated distributions
    thisyr <- yr-rings[i]
    prob <- approx(xseq, cal[i,], thisyr)$y
    segments(thisyr, i, thisyr, i+ex1*prob, col="darkgreen")
    points(thisyr, i+ex1*prob, col="darkgreen", pch=19, cex=.5)
  }

  # draw the 'wiggle match' age distribution
  these <- which(xseq <= yr)
  pol3 <- cbind(c(min(xseq[these]), xseq[these], max(xseq[these])), 
    ex2*c(0, wmd.probs[these], 0))
  polygon(pol3, col=rgb(1,0,0,.5), border=1, lwd=0.5)
  yrprob <- approx(xseq, wmd.probs, yr)$y
  points(yr, 0+ex2*yrprob, col=2, pch=19, cex=.5)
  segments(yr, 0, yr, 0+ex2*yrprob, col=2)
  
  # draw the calibration curve... 
  op <- par(new=TRUE)
  plot(0, type="n", xlim=range(xseq), xlab="", ylab="", ylim=c(920, 1350), xaxt="n", xaxs="i", yaxs="i")
  ytick <- seq(500, 1400, by=10)
  axis(side=2, at=ytick, labels = FALSE, lwd = 0.5, tck = -0.01)
  mtext("Radiocarbon Age (BP)", 2, 3, las=0)
  polygon(cc.pol, col=rgb(0,0,.8,.25), border=NA)
  lines(cc[,1], cc[,2], col=rgb(0,0,1,.25))
  axis(1, 20*(0:1e3), labels=FALSE, tck=-0.01, lwd=0.5)

  # ... and the radiocarbon dates...
  yrs <- yr-rings
  points(yrs, y, pch=19, cex=.6)
  segments(yrs, y-er, yrs, y+er)

  for(i in 1:n) {
    # ... and age distributions...
    thisyr <- yr-rings[i]
    cc.er <- approx(cc[,1], cc[,3], thisyr)$y
    uncal <- dnorm(yseq, y[i], sqrt(er[i]^2 + cc.er^2))
    ypol <- cbind(xlimplot[1]+ex3*c(0, uncal, 0), c(min(yseq), yseq, max(yseq)))
    polygon(ypol, col=rgb(0,.5,0,.25), border=rgb(0,.5,0,.25))
  }
  
    # ... and the lines (separately from above so they don't get covered by alpha green
  for(i in 1:n) {  
    thisyr <- yr-rings[i]
    mu <- approx(cc[,1], cc[,2], thisyr)$y
    hgh <- xlimplot[1] + ex3 * dnorm(approx(cc[,1], cc[,2], thisyr)$y, y[i], sqrt(er[i]^2 + cc.er^2))
    segments(thisyr, 0, thisyr, mu, col=1, lty=3, lwd=.5)
    segments(hgh, mu, -55e3, mu, col=1, lty=3, lwd=.5)
    segments(xlimplot[1], mu, hgh, mu, col=1, lwd=.5)
    segments(thisyr, mu, -55e3, mu, col=1, lty=3, lwd=.5)
    points(hgh, mu, pch=19, cex=.5, col=1)
  }
}

theseyears <- sort(unique(c(seq(min(xseq), 760, length=1e3), seq(760, 790, length=800), seq(790, max(xseq), length=1200))), decreasing=FALSE)

# requires a folder called SI_Video2, where the pngs will be placed
# this loop will take a while to run
for(i in 1:length(theseyears)) {
    png(paste0("SI_Video2/img_", 1e5+i, ".png"), height=14, width=20, res=300, units="cm")
    plot.wmd(theseyears[i])	
    dev.off()
}

# requires ffmpeg to be installed, see https://ffmpeg.org/download.html
# this command will take a while to run
system("ffmpeg -framerate 80 -pattern_type glob -y -i 'SI_Video2/*.png' -c:v libx264 -strict -2 -preset slow -pix_fmt yuv420p -vf 'scale=trunc(iw/2)*2:trunc(ih/2)*2' -f mp4 SI_Video2.mp4")
