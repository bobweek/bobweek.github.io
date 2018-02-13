n01 <- read.csv("~/Documents/Website/bobweek.github.io/2001.csv")
n16 <- read.csv("~/Documents/Website/bobweek.github.io/2016.csv")
y1 <- sum(as.numeric(n01$X.1))
y2 <- sum(as.numeric(n16$X.1))
t1 <- 2001
t2 <- 2016

lnc <- (t1*log(y2)-t2*log(y1))/(t1-t2)
a <- (log(y1)-log(y2))/(t1-t2)

c <- exp(lnc)
ts <- seq(2001,2016,by=0.1)
plot(ts,c*exp(a*ts))


m <- (y1-y2)/(t1-t2)
b <- (t1*y2-t2*y1)/(t1-t2)
plot(ts,m*ts+b)

D <- read.table("~/Documents/Website/bobweek.github.io/together_exact.txt")
plot(D$V1,D$V2,xlab="Year",ylab="Publication Count",main="Publications with keywords\n 'eco-evolution*' or 'evolutionary ecology'")
