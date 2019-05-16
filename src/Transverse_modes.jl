using Plots
plotly()

function Φ(z, λ, w0)
    zR = pi*w0^2/λ
    return atan(z/zR)
end

function w(z, λ, w0)
    zR = pi*w0^2/λ
    return w0*sqrt(1+(z/zR)^2)
end

function R(z, λ, w0)
    zR = pi*w0^2/λ
    return z+zR^2/z
end

function q(z, λ, w0)
    zR = pi*w0^2/λ
    return z + zR*1im
end

function ψ(z, p, m, λ, w0)
    zR = pi*w0^2/λ
    N = abs(m) + 2p
    return (N+1)*atan(z/zR)
end

δ(k::Integer,j::Integer) = k == j ? 1 : 0

function gLaguerre(x, p, m)
    if p == 0
        return 1
    elseif p == 1
        return 1+m-x
    else
        k = p-1
        return ((2k+1+m-x)*gLaguerre(x, k, m)-(k+m)*gLaguerre(x, k-1, m))/(k+1)
    end
end

function H(x,n)
    if n==0
        return 1
    elseif n==1
        return 2x
    else
        return 2x*H(x, n-1)-2*(n-1)*H(x, n-2)
    end
end

function r_and_θ(x, y)
    r = sqrt.(x.^2 .+ y'.^2)
    θ = atan.(y', x)
    return r, θ
end


#Laguerre-Gaussian modes
function LG(r, θ, z, p, m; λ=800e-9, w0=2e-2)
    k = 2pi/λ
    _1 = @. sqrt(2*factorial(p)/(pi*factorial(abs(m)+p)))
    _2 = @. w0/w(z, λ, w0)*(sqrt(2)*r/w(z, λ, w0))^abs(m)
    _3 = @. exp(-r^2/w(z, λ, w0)^2)
    _4 = @. gLaguerre(2r^2/w(z, λ, w0)^2,p,abs(m))
    _5 = @. exp(-1im*k*r^2/(2*R(z, λ, w0))-1im*θ*m)
    _6 = @. exp(1im*ψ(z, p, m, λ, w0))
    return  @. _1*_2*_3*_4*_5*_6
    # return _3*_4
end

#Hermite-Gaussian modes
function HG(x, y, z, l, m; λ=800e-9, w0=2e-2)
    k = 2pi/λ
    r, θ = r_and_θ(x, y)
    # outer product
    _1 =  w0/w(z, λ, w0).*H.(sqrt(2)*x/w(z, λ, w0), l).*H.(sqrt(2)*y/w(z, λ, w0), m)'
    _2 = @. exp(-(r^2)/w(z, λ, w0)^2)
    _3 = @. exp(-1im*k*(r^2)/(2*R(z, λ, w0)))
    _4 = @. exp(1im*ψ(z, l, m, λ, w0))
    return @. _1*_2*_3*_4
end



# x = LinRange{Float64}(-2,10, 100)
x = Vector(LinRange{Float64}(-5e-2,5e-2, 500))
y = Vector(LinRange{Float64}(-5e-2,5e-2, 500))
r, θ = r_and_θ(x, y)

# heatmap(x,y,real.(LG(r,θ,0,1,3)).^2)
# plot!(x,gLaguerre.(x, 2, 2))

# x = LinRange{Float64}(-2,3, 100)
# plot!(x, H.(x,4))
heatmap(x,y, real.(HG(x, y, 0, 2, 2)).^2)

# heatmap(r)
