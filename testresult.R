rm(list=ls())

###########################################3
x <- seq(0,0.15,0.03)

y50.01.F <- c(0.025, 0.1, 0.19, 0.28, 0.57, 0.725)
y50.01.M <- c(0.015, 0.05, 0.09, 0.11, 0.15, 0.28)

y100.01.F <- c(0.01, 0.09, 0.225, 0.435, 0.71, 0.9)
y100.01.M <- c(0.01, 0.03, 0.07, 0.135, 0.34, 0.535)

y50.05.F <- c(0.07, 0.16, 0.275, 0.37, 0.7, 0.85)
y50.05.M <- c(0.065, 0.145, 0.175, 0.21, 0.32, 0.475)

y100.05.F <- c(0.05, 0.2, 0.385, 0.575, 0.825, 0.94)
y100.05.M <- c(0.045, 0.12, 0.215, 0.275, 0.525, 0.72)



par(mfcol=c(2,2))

plot(x, y50.01.M , col=3, lwd=3,  xlab="c", ylab="power", ylim=c(0, 1), main=expression(paste("n=50 and ", alpha,"=0.01")), type="l")
points(x, y50.01.M , col=3, pch=16)
lines(x, y50.01.F , col=2, lwd=3)
points(x, y50.01.F , col=2, pch=16)
legend(0, 1, lwd=2, seg.len=1, legend=c("FNMEM", "NMEM"), col=c(2, 3), lty=1)

plot(x, y100.01.M , col=3, lwd=3,  xlab="c", ylab="power", ylim=c(0, 1), main=expression(paste("n=100 and ", alpha,"=0.01")), type="l")
points(x, y100.01.M , col=3, pch=16)
lines(x, y100.01.F , col=2, lwd=3)
points(x, y100.01.F , col=2, pch=16)
legend(0, 1, lwd=2, seg.len=1, legend=c("FNMEM", "NMEM"), col=c(2, 3), lty=1)

plot(x, y50.05.M , col=3, lwd=3,  xlab="c", ylab="power", ylim=c(0, 1), main=expression(paste("n=50 and ", alpha,"=0.05")), type="l")
points(x, y50.05.M , col=3, pch=16)
lines(x, y50.05.F , col=2, lwd=3)
points(x, y50.05.F , col=2, pch=16)
legend(0, 1, lwd=2, seg.len=1, legend=c("FNMEM", "NMEM"), col=c(2, 3), lty=1)

plot(x, y100.05.M , col=3, lwd=3,  xlab="c", ylab="power", ylim=c(0, 1), main=expression(paste("n=100 and ", alpha,"=0.05")), type="l")
points(x, y100.05.M , col=3, pch=16)
lines(x, y100.05.F , col=2, lwd=3)
points(x, y100.05.F , col=2, pch=16)
legend(0, 1, lwd=2, seg.len=1, legend=c("FNMEM", "NMEM"), col=c(2, 3), lty=1)