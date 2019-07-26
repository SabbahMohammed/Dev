struct Grid{T}
    τ::Array{T,1}
    dτ::T
    W::Array{T,1}
    dW::T
    ω0::T
    λ0::T
    convt::T
    convw::T
    npts::Int64
    wnpts::Int64
end

function grid(λ0, τ)
    ω0 = 2pi*c/λ0
    dτ = τ[2] - τ[1]
    dW = 2pi/(τ[end]-τ[1])
    W = Vector((0:length(τ)/2)*dW)
    convw = 1/dτ
    convt = dτ
    npts = length(τ)
    wnpts = length(W)
    return Grid(τ, dτ, W, dW, ω0, λ0, convt, convw, npts, wnpts)
end
