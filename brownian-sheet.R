
#
# build approximations of Brownian Motion and the Brownian Sheet
#

# we'll need this package to speed things up
library(JuliaCall)
julia <- julia_setup(JULIA_HOME = "/home/bb/Julia/JuliaPro-1.3.1-1/Julia/bin/")
julia_library("Distributions")
julia_library("DataFrames")

# approximate Brownian Motion using the Karhunen-Loeve Expansion on [0,1]

## r is slow
# b <- function(n,t){
#   Z = rnorm(n)
#   val = rep(0,length(t))
#   for(k in 1:n){
#     val = val + Z[k]*sin(pi*(k-0.5)*t)/((k-0.5)*pi)
#   }
#   val = sqrt(2)*val
#   return(val)
# }

## julia is fast
julia_command("
function b(n,t)
    Z = rand(Normal(),Int64(n))
    val = fill(0,length(t))
    
    for k in 1:Int64(n)
      val = val .+ Z[k].*sin.(pi.*(k-0.5).*t)./(k-0.5)
    end
    
    val = sqrt(2.0/pi).*val
    return val
end")

# test the julia function
julia_call("b", 10, 0.5)

# wrap it in R
b <- function(n,t){
  julia_call("b", n, t)
}

# have look at the approximated Brownian Motion
time <- seq(0,1,0.001)
par(bg="black",fg="white")
plot(time,b(1000,time), type='l',col="white", lwd=1,xlab="",ylab="")


# approximate Brownian Sheet with an orthogonal expansion on [0,1]

## R is slow
# W <- function(N,n,spacetime){
#   t <- spacetime$y
#   x <- spacetime$x
#   K <- length(t)
#   L <- length(x)
#   bs <- b(n,t)
#   for(i in 2:L){
#     bs <- rbind(bs,b(n,t))
#   }
#   val = matrix(0,nc=K,nr=L)
#   for(k in 0:(N-1)){
#     for(i in 1:L){
#       for(j in 1:K){
#         val[i,j] = val[i,j] + bs[i,j]*sqrt(2/pi)*sin(pi*(k+0.5)*x)/((k+0.5)*pi)
#       }
#     }
#   }
#   return(val)
# }


## julia is fast
julia_command("
function W(n,N,spacetime)
    
    t = spacetime.t
    x = spacetime.x
    K = length(t)
    L = length(x)
    
    val = zeros(K,L)
    
    for k in 1:N
      bm = b(n,t)
      for i in 1:L
        for j in 1:K
          val[i,j] += bm[j] * sin( pi * (k+0.5) * x[i] ) / (k+0.5)
        end
      end
    end
    
    val = sqrt(2.0/pi).*val
    return val
    
end")


# test by evaluating at a "point"
st <- data.frame(x=0.5,t=0.5)
julia_call("W", n, N, st)

# wrap it in R
W <- function( n, N, st ){
  julia_call("W", n, N, st)
}


################################################
#                                              #
# using plot3D to visualize the Brownian Sheet #
#                                              #
################################################

# build mesh
x <- t <- seq(0, 1, by = 0.001)
m <- data.frame(x=x,t=t)

# construct the Brownian Sheet (takes a moment)
n <- 100
N <- 100
BS <- W(n,N,m)

# plot the Brownian Sheet (also takes a moment)
require(plot3D)
require("plot3Drgl")
persp3D(z = BS, xlab = "space", bty = "bl2",
        ylab = "time", zlab = "value", clab = "depth, m",
        expand = 0.25, d = 1, phi = 20, theta = 320, resfac = 1,
        contour = list(col = "grey"),
        colkey = list(side = 1, length = 0.5))
plotrgl(lighting=TRUE)
