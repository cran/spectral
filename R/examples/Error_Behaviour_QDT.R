#### Numerischer Test für Fehler bei Quadratur demodulation

w <- 2*pi*1.125
v <- 0.1


dx <- 1e-2

x <- seq(0,1,by = dx)
# x <- sample(x = x,size = 1e2,replace = F)
N <- length(x)

x <- x[-N]
N <- length(x)

nO <- max(x)/(pi)*4*w
nO <- ifelse(round(nO) > floor(nO) ,round(nO),floor(nO))


dN <- as.integer((max(x) - nO*pi/(4*w))/dx) + 1

# windows()

res <- NULL

e <- qnorm(1-1e-3) * 2 * v / (sqrt(length(x)))


n <- 0*rnorm(length(x),0,v)
p <- 0*rnorm(1)

y <- cos(w*x+p)  + n

ak <- 2/N * sum(y*cos(w*x))
bk <- 2/N * sum(y*sin(w*x))

A <- sqrt(ak^2+bk^2)

eak <- dN/N * (1   - 2/3 * (w*max(x))^2 * (dN/N)^2 ) - 0*ifelse(ak>bk,0.5*bk/ak*(1-(1-dN/N)^2),0)


ebk <- - dN/N*( 1 - 4/3 * (w*max(x))^2 * dN^2/(N^2)) - 0*ifelse(bk>ak,0.5*ak/bk*(1-(1-dN/N)^2),0)
# das geht aber im Quotienten dürfte kein n stehen.

e <- sqrt(ebk^2+eak^2)

A2 <- sqrt((ak/(1 + eak))^2 + (bk/(1+ebk))^2)

plot(x,y)
curve(cos(w*x+p),add=T,n = length(x))
curve(cos(2*w*x+p), from = 0,to=max(x),add=T)
abline(h=0,v=nO*2*pi/w,lty=3)
