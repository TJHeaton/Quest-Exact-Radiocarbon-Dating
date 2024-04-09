# load an R package which has the IntCal20 calibration curve
if(!'rintcal'%in% installed.packages())
  install.packages('rintcal')
require(rintcal) 
cc <- ccurve()
cc[,4] <- cc[,1] - 1950 # work in BC (AD negative)

xlim <- c(1880, 1390) # manual, BC
ylim <- c(3000, 3600) # manual
xseq <- seq(xlim[1], xlim[2], length=1000)
xseq <- sort(unique(xseq, seq(1760, 1520, length=3e3)), decreasing=TRUE) # zoom in
yseq <- seq(ylim[1], ylim[2], length=length(xseq))

y <- 3350
er <- 10

# the raw C14 age, likelihood assuming a normal distribution
l.yseq <- dnorm(yseq, y, er)

# calibration
mu.xseq <- approx(cc[,4], cc[,2], xseq)$y # BC
sd.xseq <- approx(cc[,4], cc[,3], xseq)$y
l.xseq <- dnorm(y, mu.xseq, sqrt(er^2 + sd.xseq^2)) # calibration
l.xseq <- l.xseq/sum(l.xseq)

cc.sel <- cbind(xseq, mu.xseq, sd.xseq) # interpolated calibration curve, in BC

# exaggeration factors for the distributions
ex1 <- .3/max(l.xseq)*(ylim[2] - ylim[1]) # calibrated distribution
ex2 <- 1500 # does not depend on other values

# now place on plot coordinates
xax.offset <- 2
yax.offset <- 4
xpol <- cbind(c(xlim[1], xseq, xlim[2]), ylim[1] + ex1*c(0, l.xseq, 0))
ypol <- cbind(xlim[1] - ex2*c(0, l.yseq, 0), c(ylim[1], yseq, ylim[2]))
ccpol <- cbind(c(xseq, rev(xseq)), c(mu.xseq-sd.xseq, rev(mu.xseq+sd.xseq)))

# as function of xseq[i]
plot.cal <- function(i) {
  op <- par(xaxs="i", yaxs="i", mar = c(5, 4.8, 2, 2) + 0.1)
  plot(cc.sel[,1], cc.sel[,2], type="l", xlim=xlim, ylim=ylim,
    xlab="Calendar Age (cal BC)", ylab="", cex.lab=1, las=1, col=rgb(213,94,0, maxColorValue = 255))
  mtext("Radiocarbon Age (BP)", side = 2, line = 3.2, cex = 1.3)
  legend("topright", legend = "IntCal20", 
    lty = c(-1), pch = c(15), lwd = 2,
    col = rgb(213,94,0, alpha = 127, maxColorValue = 255), cex = 1, pt.cex = 2)
  xtick <- seq(xlim[2], xlim[1], by = 10)
  axis(side = 1, at = xtick, labels = FALSE, lwd = 0.5, tck = -0.015)
  ytick <- seq(ylim[1], ylim[2], by = 10)
  axis(side = 2, at = ytick, labels = FALSE, lwd = 0.5, tck = -0.015)
  
  # draw the calibration curve
  polygon(ccpol, col=rgb(213,94,0, alpha = 127, maxColorValue = 255), border=NA)

  # use y.pol, i.e., static, not adapting to changing cc.sd
  polygon(ypol, col=rgb(1,0,1,.5))    

  # evolving probability on radiocarbon scale
  mu.x <- approx(cc[,4], cc[,2], xseq[i])$y
  l.y <- dnorm(mu.x, y, er)
  points(xlim[1] - ex2 * l.y, mu.xseq[i], pch=19, cex=.8)
  
  segments(99e3, mu.xseq[i], xlim[1] - ex2 * l.xseq[i], mu.xseq[i])
 
  # calibrated distribution
  polygon(rbind(xpol[1:i,], c(xpol[i], ylim[1])), col=rgb(1,0,1,.2))
  points(xpol[i], ylim[1] + ex1 * l.xseq[i], pch=19, cex=.8)
  segments(xpol[i], ylim[1], xpol[i], ylim[1] + ex1 * l.xseq[i], lwd=2)

  # connecting lines
  segments(xseq[i], -99e3, xpol[i], mu.xseq[i], lty=2, col=rgb(0,0,0,.5))
  segments(99e3, mu.xseq[i], xseq[i], mu.xseq[i], lty=2, col=rgb(0,0,0,.5))
  
  # draw the radiocarbon date
  xpos <- xlim[1]-2*xax.offset
  points(xpos, y, pch=20, cex=.5)
  segments(xpos, y-er, xpos, y+er, lwd=2)
}

scalefac <- 1.5
# requires a folder called SI_Video1, where the pngs will be placed
# this loop will take a while to run
for(i in 1:length(xseq)) {
  png(paste0("SI_Video1/img_", 2e5+i, ".png"), width = scalefac * 7, height = scalefac * 3, units = 'in', res = 300)
  plot.cal(i)
  dev.off()
}

# requires ffmpeg to be installed, see https://ffmpeg.org/download.html
# this command will take a while to run
system("ffmpeg -framerate 40 -pattern_type glob -y -i 'SI_Video1/*.png' -c:v libx264 -strict -2 -preset slow -pix_fmt yuv420p -aspect 7:3 -vf 'scale=2560:-2' -f mp4 SI_Video1.mp4")
