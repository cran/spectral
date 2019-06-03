# create two sin-functions
x_orig <- seq(0,1,by=1e-2)
y_orig <- sin(10*2*pi*x_orig) + 0.1*sin(2*2*pi*x_orig)

# make a 10% gap
i <- round(length(x_orig)*0.2) : round(length(x_orig)*0.3)
x <- x_orig
y <- y_orig
x[i] <- NA
y[i] <- NA


# calculating the lomb periodogram
l <- spec.lomb(x = x, y = y,ofac = 20,mode = "normal")

# select a frequency range
m <- rbind(c(9,11))
# select and reconstruct the most significant component
l2 = filter.lomb(l, x_orig, filt = m)

# plot everything
par(mfrow=c(2,1),mar = c(4,4,2,4))
plot(x,y,"l", main = "Gapped signal")
lines(l2$x, l2$y,lty=2)
legend("bottomleft",c("gapped","10Hz component"),lty=c(1,2))

plot(l,main = "Spectrum")


### Multivariate -- 3D Expample ###
require(lattice)
fx <- 8.1
fy <- 5
fz <- 2

# creating frequency space
f <- expand.grid( fx = seq(-10,10,by = 0.25)
                  ,fy = seq(-10,10,by = 0.25)
                  ,fz = 0:3
)

# creating spatial space
pts <- expand.grid( x = seq(0,1,by = 0.025)
                   ,y = seq(0,1,by = 0.025)
                   ,z = seq(0,1,by = 0.025)
)

# gapping 30%
i <- sample(1:dim(pts)[1],0.7*dim(pts)[1])
pts <- pts[i,]

# caluculating function
pts$val <- cos(2*pi*(  fx*pts$x
                     + fy*pts$y
                     + fz*pts$z
                    ) + pi/4
              )

# display with lattice
levelplot(val~x+y,pts,subset = z == 0,main = "with z = 0")

# calculating lomb takes a while
# or we sample only a few points
# which enlarges the noise but accelerates the calculation
l <- spec.lomb(y = pts[sample(1:dim(pts)[1],5e2),]
               ,f = f
               ,mode = "generalized"
               )

# name the stripes
l$fz_lev <- factor(x = paste("fz =",l$fz)
)

# display output
levelplot(PSD~fx+fy|fz_lev,l)

