function gnlse(pulse::Pulse, fiber::Fiber, betas::Array{Any,1}; dispersion::Bool = true,  fr::Float64=0.0, shock::Bool=false, plasma::Bool=false, nsaves::Int64=1000)
    ω0 = 2pi*c/pulse.λ0
    γ = getGamma(fiber.gas, fiber.p, fiber.Tk, ω0, fiber.diameter)

    # grid
    τ = pulse.τ
    dτ = τ[2] - τ[1]
    convt = length(τ)*dτ
    W = 2pi.*fftshift(fftfreq(length(τ), 1/dτ))
    dW = W[2] - W[1]
    convw = dW/(2pi)

    # define the indices for where the signal after them ignored
    # the limits here are 200 nm - 2000 nm
    # downlim = findfirst(x -> x>(2pi*3e8/downlim - ω0), W)
    # uplim = findfirst(x -> x>(2pi*3e8/uplim - ω0), W)

    #dispersion
    if dispersion == true
        if length(betas) > 0

            β2 = wbeta2(ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2)
            betas0 = [β2]
            β = zeros(Float64, size(W))
            for i in 1:length(betas0)
                @. β = β + betas0[i]/factorial(i+1)*(W)^(i+1)
            end
            L = 1im*β
        else
            L = @. 1im*L_operator.(W,ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - fiber.α/2
        end

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
    τ1 = 62.5e-15
    τ2 = 77e-15
    RW = sdo(τ, τ1, τ2)

    #windows
    gw(z) = (convt) .* exp.(-0.5 .*((W.+ 2e15 .-  maximum(W)./2)./(maximum(W).*0.5)).^30).* exp.(-z .* L)
    gt(z) =  exp.(-0.5 .*((τ)./(maximum(τ).*0.5)).^30) .* convw

    Uw0 = shifted_ifft(pulse.data) .* convt
    zs = (0.0, fiber.length)
    Z = LinRange(0, fiber.length, nsaves)
    prob = ODEProblem(rhs, Uw0, zs, (L, Progress(floor(Int, Float64(fiber.length*1000)), 2),
                      W, τ, ω0, γ, RW, fiber, pulse, fr, shock, plasma, dW, dτ, convw, convt, gw, gt))
    sol = solve(prob, Vern7(), reltol=1e-4, abstol=1e-14, saveat=fiber.length/nsaves)
    res = zeros(Complex{Float64}, size(sol))
    @simd for i in eachindex(sol.t)
        @inbounds res[:,i] = Complex{Float64}.(sol[:,i] .* exp.(L .* sol.t[i]))
    end
    rest = zeros(Complex{Float64}, size(sol))
    @simd for i in eachindex(sol.t)
        @inbounds rest[:,i] = Complex{Float64}.(shifted_fft(res[:,i])).*convw
    end

    τ, W, rest, res, Z
end

function rhs(dUw, Uw, para, z)
    L = para[1]
    W, τ, ω0, γ, RW, fiber, pulse, fr, shock, plasma, dW, dτ, convw, convt, gw, gt = para[3:end]

    update!(para[2], max(floor(Int, Float64(z*1000)), 1))
    U = shifted_fft(Uw .* exp.(L .* z)) .* convw
    Uw = shifted_ifft(abs2.(U) .* U).*convt
    if fr != 0.0
        covRaman(U, RW, convt, convw, fr)
    end
    if shock == true
        @. Uw *= (1.0 .+ W./ω0)
    end

    if plasma == true
        Et = getElectricField(pulse, fiber)
        #Pl = plasmaTerms(Et, τ, U, fiber.Ip, fiber.p, pulse.λ0, fiber.diameter).* (length(τ)*dτ)
        Pl = geisslerPlasma(Et, τ, U, fiber.Ip, fiber.p, pulse.λ0, fiber.diameter, convw, convt).* (length(τ)*dτ)
        #Pl = fedotovPlasma(Et, τ, U, fiber.Ip, fiber.p, pulse.λ0, fiber.diameter).* (length(τ)*dτ)
        # @. dUw = (1im*γ*uuw + Pl) * exp(-z * L)
        dU =  @. (1im*γ*Uw + Pl)
        dUt = shifted_fft(dU) .* gt(z)
        dUw .= shifted_ifft(dUt).* gw(z)
    else
        # @. dUw = (1im*γ*uuw) * exp(-z * L)
        dU = @. (1im*γ*Uw)
        dUt = shifted_fft(dU) .* gt(z)
        dUw .= shifted_ifft(dUt).* gw(z)
    end

end

#based on soliton order
function gnlse(order::Float64, τ::Array{Float64,1}, betas::Array{Float64,1}, λ0::Float64, FWHM::Float64, fiber::Fiber; dispersion::Bool = true,  fr::Float64=0.0, shock::Bool=false, plasma::Bool=false, nsaves=1000)
    τ0 = FWHM/1.76
    ω0 = 2pi*c/λ0
    γ = getGamma(fiber.gas, fiber.p, fiber.Tk, ω0, fiber.diameter)

    P0 = abs(betas[1])*(order^2)/(γ*τ0^2)
    energy = 2τ0*P0
    pulse = sechPulse(τ, energy, FWHM, λ0)

    # grid
    τ = pulse.τ
    dτ = τ[2] - τ[1]
    convt = length(τ)*dτ
    W = 2pi.*fftshift(fftfreq(length(τ), 1/dτ))
    dW = W[2] - W[1]
    convw = (dW/(2pi))



    #dispersion
    if dispersion == true
        β = zeros(Float64, size(W))
        for i in 1:length(betas)
            @. β = β + betas[i]/factorial(i+1)*(W)^(i+1)
        end
        L = 1im*β
    else
        L = 0
    end

    #Raman
    if fr != Nothing
        γ = γ/(1-fr)
        Ht = Htot(τ)
        Ht = f2(Ht, τ)
        RT = Ht
        RT = RT/trapz(τ, RT)
        RW = shifted_ifft(RT)*convt
    end
    function rhs(dUw, Uw, para, z)
        L = para[1]
        # println(z)
        update!(para[2], max(floor(Int, Float64(z*1000)), 1))
        U = shifted_fft(Uw .* exp.(L .* z)) .* convw
        uu = abs2.(U) .* U
        uuw = shifted_ifft(uu).*convt
        if fr != Nothing
            # println(fr)
            Iw = shifted_ifft(abs2.(U)).*convt
            M = shifted_fft(Iw.*RW) .* U .* convw
            Mw = shifted_ifft(M) .* convt
            @. uuw = (1.0 - fr) .* uuw .+ fr.* Mw
        end
        if shock == true
            @. uuw = uuw + uuw * W/ω0
        end

        if plasma == true
            Et = getElectricField(pulse, fiber)
            Pl = plasmaTerms(Et, τ, U, fiber.Ip, fiber.p, pulse.λ0, fiber.diameter).* (length(τ)*dτ)

            # savefig(plot(τ, abs2.(plasma)), "plas")
            @. dUw =  (1im*γ*uuw + Pl) * exp(-z * L)
        else
            @. dUw = (1im*γ*uuw) * exp(-z * L)
        end

    end
    Uw0 = shifted_ifft(pulse.data) .* convt
    zs = (0.0, fiber.length)

    prob = ODEProblem(rhs, Uw0, zs, (L, Progress(floor(Int, Float64(fiber.length*1000)), 2)))
    sol = solve(prob, Vern7(), reltol=1e-8, abstol=1e-16, saveat=fiber.length/nsaves)
    res = zeros(Complex{Float64}, size(sol))
    for i in eachindex(sol.t)
        res[:,i] = Complex{Float64}.(sol[:,i] .* exp.(L .* sol.t[i]))
    end
    rest = zeros(Complex{Float64}, size(sol))
    for i in eachindex(sol.t)
        rest[:,i] = Complex{Float64}.(shifted_fft(res[:,i])).*convw
    end

    τ, W.+ω0, rest, res
end

function myplot(τ, W, solt, sol, name)
    ldw = 10.0.*log10.(abs2.(sol) .+ 1e-27);

    ldt = abs2.(solt);

    savefig(heatmap(ldw, clim=(maximum(ldw)-50.0, maximum(ldw)), color = :pu_or), "$name_sol")

    savefig(heatmap(ldt, clim=(maximum(ldt)-70.0, maximum(ldt)), color = :pu_or), "$name_solt")
end
