function uppe(pulse::Pulse, fiber::Fiber, betas::Array{Float64,1}; dispersion::Bool = true, thg::Bool = false,  fr::Float64=0.0, shock::Bool=false, plasma::Bool=false, nsaves::Int64=1000)
    ω0 = 2pi*c/pulse.λ0
    # grid
    τ = pulse.τ
    dτ = τ[2] - τ[1]
    convt = length(τ)*dτ
    dW = 2pi/(τ[end]-τ[1])
    W = (0:length(τ)/2)*dW
    convw = 1/dτ
    convt = dτ
    #convw = dW/(2pi)
    # if dispersion == true
    #     L = @. 1im*L_operator.(W,ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - fiber.α/2
    #     # L = []
    #     # arr = readdlm("L_operator.txt")
    #     # for i in 1:length(arr)
    #     #     append!(L, parse(Complex{Float64}, arr[i]))
    #     # end
    # else


    #dispersion
    if dispersion == true
        L = @. 1im*LE_operator.(W,ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - fiber.α/2
        # L = []
        # arr = readdlm("L_operator.txt")
        # for i in 1:length(arr)
        #     append!(L, parse(Complex{Float64}, arr[i]))
        # end
    else
    # if dispersion == true
    #     β = zeros(Float64, size(W))
    #     for i in 1:length(betas)
    #         @. β = β + betas[i]/factorial(i+1)*(W)^(i+1)
    #     end
    #     L = 1im*β
    # else
        L = 0
    end
    #Raman
    #RW = raman(fr, γ, convt)

    #windows
    gw = exp.(-100 .*((W .-  maximum(W)./2.2)./(maximum(W).*0.47)).^40)
    gt =  exp.(-0.5 .*((τ)./(maximum(τ).*0.8)).^30)

    #some calculation done insde rhs
    wb0 = W.^2.0./(2*c^2.0.*ϵ0.*wbeta0.(W,fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2))
    wbdiff = wb0[end] - wb0[1]
    wdiff = W[end] - W[1]
    wb = wbdiff/wdiff*W

    # pulse2 = gaussianPulse(τ, 0.1e-6, 30e-15, 400e-9);
    Et = getElectricFieldE(pulse, fiber)# .+ getElectricFieldE(pulse2, fiber)
    Ew0 = rfft(Et)*dτ

    #pre-allocate arrays in the memory
    ret = zeros(Float64, size(Et))#
    Pe0 = zeros(Float64, size(Et))
    Pe = zeros(Float64, size(Et))
    Na = zeros(Float64, size(Et))
    Ne = zeros(Float64, size(Et))
    WE = zeros(Float64, size(Et))
    Pe_l = zeros(Float64, size(Et))
    Pe_f = zeros(Float64, size(Et))
    PNL = zeros(Complex{Float64}, size(Ew0))
    dEt = zeros(Float64, size(Et))

    zs = (0.0, fiber.length)
    Z = LinRange(0, fiber.length, nsaves)
    prob1 = ODEProblem(rhs2, Ew0, zs, (L, Progress(floor(Int, Float64(fiber.length*1000)), 2),
                      W, τ, ω0, fiber, pulse, gw, gt, convw, convt, Et, ret, PNL, dEt, Pe,
                      Pe0, Na, Ne, WE, Pe_l, Pe_f, wb, plasma, thg))
    # sol = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-14, force_dtmin=true, saveat=fiber.length/nsaves)
    sol = solve(prob1, Tsit5(), reltol=1e-4, abstol=1e-10, saveat=fiber.length/nsaves)
    res = zeros(Complex{Float64}, size(sol))
    @simd for i in eachindex(sol.t)
        @inbounds res[:,i] = Complex{Float64}.(sol[:,i] .* exp.(-L .* sol.t[i]))#
    end
    # rest = zeros(Float64, (length(τ), nsaves+1))
    # @simd for i in eachindex(Z)
    #     @inbounds rest[:,i] = Float64.(irfft(res[:,i], length(τ)) .* 1/dτ)
    # end
    # rest = zeros(Complex{Float64}, (length(τ), nsaves+1))
    # @simd for i in eachindex(Z)
    #     Ew = 2.0.*sol[:,i]
    #     Ew = append!(zeros(Complex{Float64}, (length(τ)- length(sol[:,i]))), sol[:,i])
    #     @inbounds rest[:,i] = Complex{Float64}.(ifft(Ew)).*1/dτ
    # end
    rest = zeros(Complex{Float64}, (length(τ), nsaves+1))
    @simd for i in eachindex(Z)
        Ew = res[:,i].*2
        Ew1 = append!(zeros(Complex{Float64}, (length(τ)- length(Ew))), Ew)
        @inbounds rest[:,i] = Complex{Float64}.(ifft(Ew1))
    end

    τ, W, rest, res, Z, L
end

function rhs2(dEw, Ew, para, z)
    L = para[1]
    W, τ, ω0, fiber, pulse, gw, gt, convw, convt, Et, ret, PNL, dEt,
        Pe, Pe0, Na, Ne, WE, Pe_l, Pe_f, wb, plasma, thg = para[3:end]
    dτ = τ[2] - τ[1]
    dW = 2pi/(τ[end]-τ[1])
    #print(z)
    update!(para[2], max(floor(Int, Float64(z*1000)), 1))

    #=
    why the scaling is like this?
    convw = 1/dτ
    convt = dτ
    from the fact that how DFT work (https://en.wikipedia.org/wiki/Discrete_Fourier_transform)
    it is clear that the algorithm does not include the integral step (ie dt), so we multiply fft
    by dτ and ifft by df*N which gives 1/dt(the *N is because in the algorithm it divide by N)
    =#

    Ew = Ew.*exp.(-L.*z)

    Et .= irfft(Ew, length(τ)) .* convw # real.(ifft(cat(conj(Ew[end-1:-1:2]), Ew, dims =1))) .* convw

    if plasma == true
        Pe .= geisslerPlasmaUPPE(Et, τ, dτ, fiber.Ip, fiber.p, Pe0, Na, Ne, WE, Pe_l, Pe_f)
        ret .= χ3(fiber.gas, fiber.p, fiber.Tk)*ϵ0.*Et.^3 .+ Pe
    elseif thg == true
        ret .= χ3(fiber.gas, fiber.p, fiber.Tk)*ϵ0.*(Et).^3
    else
        ret .= χ3(fiber.gas, fiber.p, fiber.Tk)*ϵ0.*abs2.(Et).*Et
    end

    PNL .= -1im.*wb.*rfft(ret) .* convt

    dEt .= irfft(PNL, length(τ)) .* convw .* gt

    # println(z)
    dEw .= rfft(dEt) .*exp.(L.*z).* convt.* gw
end
#
# function testrhs(pulse, fiber, z)
#     τ = pulse.τ
#     dτ = τ[2] - τ[1]
#     dW = 2pi/(τ[end]-τ[1])
#     ω0 = 2pi*c/pulse.λ0
#     W = (0:length(τ)/2)*dW
#
#
#
#     convw = dW/(2pi)
#     convt = length(τ)*dτ
#
#     L = @. 1im*LE_operator.(W,ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - fiber.α/2
#
#     Et = getElectricFieldE(τ, pulse.data, fiber.diameter, ω0)
#
#     #Et = getElectricFieldE(τ, pulse.data, fiber.diameter, ω0)
#     Ew = ifft(Et) .* convt
#
#     #Ew = Ew[2^13:end].*exp.(L.*z)
#
#     Et1 = real.(fft(Ew)) .* convw
#
#     prob1 = init(rhs2, Ew0, zs)
#
#     # W0 = 2pi.*fftshift(fftfreq(length(τ), 1/dτ))
#     # Uw0 = shifted_ifft(pulse.data) .* convt
#
#
#     # Ew = ifft(Et) .* convt
#
#     # println(NumericalIntegration.integrate(W0, abs2.(Ew)))
#     # println(NumericalIntegration.integrate(W0, abs2.(Uw0)))
#
#     # #Et = real.(ifft(cat(conj(Ew[end-1:-1:2]), Ew, dims =1))) .* convw
#
#     # Ew = cat(conj(Ew[end-1:-1:2]), Ew, dims =1);
#
#     #Et = real.(fft(cat(conj(Ew[end-1:-1:2]), Ew, dims =1))) .* convw
#     #Et = real.(ifft(cat(Ew0, conj(Ew0[end-1:-1:2]), dims =1))).* convw
#
#     # Ew0 = rfft(Et).* convt
#
#     # ret = cat(conj(Ew0[end-1:-1:2]), Ew0, dims =1)
#
#
#     return τ, Et, Et1
# end
#
# function uppe2(pulse::Pulse, fiber::Fiber, betas::Array{Float64,1}; dispersion::Bool = true,  fr::Float64=0.0, shock::Bool=false, plasma::Bool=false, nsaves::Int64=1000)
#     ω0 = 2pi*c/pulse.λ0
#     #γ = getGamma(fiber.gas, fiber.p, fiber.Tk, ω0, fiber.diameter)
#
#     # grid
#     τ = pulse.τ
#     dτ = τ[2] - τ[1]
#     convt = length(τ)*dτ
#     dW = 2pi/(τ[end]-τ[1])
#     W = (0:length(τ)/2)*dW
#     #convw = dW/(2pi)
#     # if dispersion == true
#     #     L = @. 1im*L_operator.(W,ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - fiber.α/2
#     #     # L = []
#     #     # arr = readdlm("L_operator.txt")
#     #     # for i in 1:length(arr)
#     #     #     append!(L, parse(Complex{Float64}, arr[i]))
#     #     # end
#     # else
#
#
#     #dispersion
#     if dispersion == true
#         L = @. 1im*LE_operator.(W,ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - fiber.α/2
#         # L = []
#         # arr = readdlm("L_operator.txt")
#         # for i in 1:length(arr)
#         #     append!(L, parse(Complex{Float64}, arr[i]))
#         # end
#     else
#     # if dispersion == true
#     #     β = zeros(Float64, size(W))
#     #     for i in 1:length(betas)
#     #         @. β = β + betas[i]/factorial(i+1)*(W)^(i+1)
#     #     end
#     #     L = 1im*β
#     # else
#         L = 0
#     end
#     #Raman
#     #RW = raman(fr, γ, convt)
#
#     #windows
#     gw(z) = exp.(-0.5 .*((W.+ 2e15 .-  maximum(W)./2)./(maximum(W).*0.5)).^30).* exp.(-z .* L)
#     gt(z) =  exp.(-0.5 .*((τ)./(maximum(τ).*0.5)).^30)
#
#     Et = getElectricFieldE(τ, pulse.data, fiber.diameter, ω0)
#     Ew0 = rfft(Et)
#     zs = (0.0, fiber.length)
#     Z = LinRange(0, fiber.length, nsaves)
#     prob1 = ODEProblem(rhs2, Ew0, zs, (L, Progress(floor(Int, Float64(fiber.length*1000)), 2),
#                       W, τ, ω0, fiber, pulse))
#     integrator = init(prob1, Vern7(), reltol=1e-4, abstol=1e-14)
#     return integrator
# end
