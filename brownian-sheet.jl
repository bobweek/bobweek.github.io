using Distributions

using AbstractPlotting, GLMakie

function b(n,t)
    Z = rand(Normal(),Int64(n))
    val = fill(0,length(t))

    for k in 1:Int64(n)
      val = val .+ Z[k].*sin.(pi.*(k-0.5).*t)./(k-0.5)
    end

    val = sqrt(2.0/pi).*val
    return val
end

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

end

BrSh = W(2000,2000,(x=0:0.001:1,t=0:0.001:1))

BrSh *= 10
scene = Scene(backgroundcolor = :black)
surface!(BrSh; shading=false, show_axis=false, colormap = :deep)
scene
