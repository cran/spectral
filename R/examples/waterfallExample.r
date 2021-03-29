#### noisy signal with amplitude modulation ####
x <- seq(0,3, length.out = 1000)
# original data
# extended example from envelope function
y <- 1*(abs(x-1.5))*sin(10*2*pi*x) + ifelse(x > 1.5,sin(15*(1+0.25*(x - 1.5))*2*pi*x),0)
ye <- base::Re(envelope(y))

par(mfrow=c(2,1),mar=c(1,3.5,3,3),mgp=c(2.5,1,0))

# plot results
plot(x,y,type="l",lwd=1,col="darkgrey",lty=2,ylab="y",main="Original Data",xaxt="n",xlab="")
lines(x,ye)
legend("bottomright",c("modulated","envelope"),col=c("grey","black"),lty=c(2,1))

par(mar=c(3.5,3.5,2,0))
wf <- waterfall(y,x,nf = 3)
# rasterImage2(x = wf$x, y = wf$fx, z = wf$A
#              ,ylim = c(0,60))

plot(wf,ylim=c(0,40),main="Waterfall")


#### uncertainty principle ####
#
# take a look at the side effects
# at [0,30] and [1,0]
#
# With a large steepness e.g. n = 50 you will gain
# artefacts.
#
# if frequency is not stationary
# PSD becomes > 1 depending on the type of band filter.
#
###############################
x <- seq(0,1, length.out=1500)
y <- sin(100*x*x)

FT <- spec.fft(x = x, y = y)
wf <- waterfall(y,x)

par(mfrow=c(2,1),mar=c(1,3.5,3,3),mgp=c(2.5,1,0))
# plot results
plot(x,y,type="l",lwd=1,col="darkgrey",lty=2,ylab="y",main="Original Data",xaxt="n",xlab="")

par(mar=c(3.5,3.5,2,0))
plot(wf
             ,ylim=c(0,40),main="Waterfall"
             )
abline(h = 25, lty = 3, lwd = 3, col = "grey")
range(wf$PSD,na.rm = TRUE)
range(wf$A)

###### effect of missing values #####
#
# 10% random missing values cause a
# distortion and a miss scaling of
# the PSD value, which becomes >1 now.
# This depends on the type of band pass
# filter selected.
#
#####################################
x <- seq(0,5, length.out=500)
y <- sin(2*pi * 15 * x + 2*1*cos(2*pi*0.5*x))

# delete 10% of the data
y[sample(length(y),size = 50)] <- NA

wf <- waterfall(y,x,type = "b")

par(mfrow=c(2,1),mar=c(1,3.5,3,3),mgp=c(2.5,1,0))
# plot results
plot(x,y,type="l",lwd=1,col="darkgrey",lty=2,ylab="y",main="Original Data",xaxt="n",xlab="")

par(mar=c(3.5,3.5,2,0))
plot(wf
     ,ylim=c(10,20),main="Waterfall"
)
abline(h = 25, lty = 3, lwd = 3, col = "grey")

# check the PSD range
range(wf$PSD)
range(wf$A)
