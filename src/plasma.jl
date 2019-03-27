function ionRateAdk(Et, Ip)
    Ip *= e
    ns = sqrt(13.59844*e/Ip)
    wp = Ip/hbar
    cn2 = 2^(2*ns)/(ns*gamma(ns+1)*gamma(ns))
    wt = @. e*abs(Et)/sqrt(2*m_e*Ip)
    WE = @. wp*cn2*(4wp/wt)^(2ns-1)*exp(-4.0/3.0*wp/wt)
    WE[notfinite.(WE)] .= 0.0
    WE
end

function adk(Et, τ, Ip, p)
    WE = ionRateAdk(Et, Ip)
    Na = exp.(-cumul_integrate(τ, WE))
    Ne = p*N0.*(1.0.-Na)
    Na = p*N0.*Na
    Na[notfinite.(Na) .| (Na .< 0.0)] .= 1.0
    Ne[notfinite.(Ne) .| (Ne .< 0.0)] .= 0.0
    [Na, Ne, WE]
end

function plasmaTerms(Et, τ, U, Ip, p, λ0, diameter)
    ω0 = 2pi*c/λ0
    dτ = τ[2] - τ[1]
    Na, Ne, WE = adk(Et, τ, Ip, p)
    frequancy_term = @. e^2.0*Ne/(2ϵ0*m_e*c*ω0)
    Aeff = effectiveArea(diameter)
    loss_term = @. Na*WE*Aeff*Ip*e/(2.0*abs2(U))
    loss_term[notfinite.(loss_term)] .= 0.0
    Pl = -1im*shifted_ifft(U .* frequancy_term)  .- shifted_ifft(U .* loss_term)
    return Pl
end

function geisslerPlasma(Et, τ, U, Ip, p, λ0, diameter, convw, convt)
    dτ = τ[2] - τ[1]
    ω0 = 2pi*c/λ0
    Aeff = effectiveArea(diameter)
    Na, Ne, WE = adk(Et, τ, Ip, p)
    #loss_term = 1im*ω0*Ip*e/(c*ϵ0*Aeff) .* cumul_integrate(τ, Na.*WE./Et)
    #loss_term = 1im*ω0*Ip*e/(c^2*ϵ0^2*Aeff) .* cumul_integrate(τ, Na.*WE./Et)
    #loss_term = 1im*ω0*Ip*e/(c*ϵ0*sqrt(Aeff)) .* cumul_integrate(τ, Na.*WE./Et)
    loss_term = 1im*ω0*Ip*e/(c^2*ϵ0*sqrt(Aeff)) .* cumul_integrate(τ, Na.*WE./Et)
    loss_term[notfinite.(loss_term)] .= 0.0
    frequancy_term = cumul_integrate(τ, Ne.*Et)
    #frequancy_term = 1im*ω0*e^2/(m_e*c^2*ϵ0^2*Aeff).*cumul_integrate(τ, frequancy_term)
    #frequancy_term = 1im*ω0*e^2/(m_e*c^2*ϵ0^2*Aeff).*cumul_integrate(τ, frequancy_term)
    frequancy_term = 1im*ω0*e^2/(m_e*c^2*ϵ0*sqrt(Aeff)).*cumul_integrate(τ, frequancy_term)

    # Hilbert transform 5.19 in notes
    loss_term = 2.0.*shifted_ifft(loss_term) .*convt
    loss_term[τ.<0] .= 0.0
    frequancy_term = 2.0.*shifted_ifft(frequancy_term) .*convt
    frequancy_term[τ.<0] .= 0.0
    loss_term = shifted_fft(loss_term) .* convw
    frequancy_term= shifted_fft(frequancy_term) .* convw

    # back to envelop 5.23
    loss_term .= exp.(1im.*ω0 .* τ)
    frequancy_term .= exp.(1im.*ω0 .* τ)

    return (shifted_ifft(frequancy_term) .+ shifted_ifft(loss_term)).*convt
end

function geisslerPlasmaUPPE(Et, τ, dτ, Ip, p, Pe0, Na, Ne, WE, Pe_l, Pe_f)
    Na, Ne, WE = adk(Et, τ, Ip, p)
    Pe_l .= Ip*e .* cumul_integrate(τ, Na.*WE./Et)
    Pe_f .= cumul_integrate(τ, Ne.*Et)
    Pe_f .= e^2/m_e.*cumul_integrate(τ, Pe_f)
    Pe0 .= Pe_f .+ Pe_l
    Pe0[notfinite.(Pe0)] .= 0.0
    return Pe0
end

# function geisslerPlasma(Et, τ, U, Ip, p, λ0, diameter)
#     dτ = τ[2] - τ[1]
#     ω0 = 2pi*c/λ0
#     Aeff = effectiveArea(diameter)
#     Na, Ne, WE = adk(Et, τ, Ip, p)
#     loss_term = @. Na*WE*Aeff*Ip*e/(2.0*abs2(U))
#     loss_term[notfinite.(loss_term)] .= 0.0
#     frequancy_term = cumul_integrate(τ, Ne.*Et)
#     frequancy_term = -1im*e^2/(2*m_e*c*ϵ0*sqrt(Aeff)).*cumul_integrate(τ, frequancy_term)
#     return  (shifted_ifft(frequancy_term) .- shifted_ifft(loss_term .* U)) # does it have to be multiplied by U ??
# end

function fedotovPlasma(Et, τ, U, Ip, p, λ0, diameter)
    ω0 = 2pi*c/λ0
    dτ = τ[2] - τ[1]
    Na, Ne, WE = adk(Et, τ, Ip, p)
    Aeff = effectiveArea(diameter)
    frequancy_term = @. 2pi*e^2.0*Ne/(c*m_e*ω0)
    loss_term = @. Na*WE*Aeff*Ip*e/(2.0*abs2(U))
    loss_term[notfinite.(loss_term)] .= 0.0
    return -1im*shifted_ifft(U .* frequancy_term)  .- shifted_ifft(U .* loss_term)
end

# function fedotovPlasma(Et, τ, U, Ip, p, λ0, diameter)
#     ω0 = 2pi*c/λ0
#     dτ = τ[2] - τ[1]
#     Na, Ne, WE = adk(Et, τ, Ip, p)
#     Aeff = effectiveArea(diameter)
#     frequancy_term = @. 2pi*e^2.0*Ne/(c*m_e*ω0*Aeff)
#     loss_term = @. Na*WE*Aeff*Ip*e/(2.0*abs2(U))
#     loss_term[notfinite.(loss_term)] .= 0.0
#     return -1im*shifted_ifft(U .* frequancy_term)  .- shifted_ifft(U .* loss_term)
# end

#*sqrt(ϵ0*c*Aeff/2)

function geisslerPlasmaE(Et, τ, Ip, p)
    #dτ = τ[2] - τ[1]
    #ω0 = 2pi*c/λ0
    #Aeff = effectiveArea(diameter)
    Na, Ne, WE = adk(Et, τ, Ip, p)
    loss_term = Ip*e.*cumul_integrate(τ, Na*WE*/(Et))
    loss_term[notfinite.(loss_term)] .= 0.0
    loss_term[τ.<0] .= 0.0
    frequancy_term = cumul_integrate(τ, Ne.*Et)
    frequancy_term = -1im*e^2/(m_e).*cumul_integrate(τ, frequancy_term)
    frequancy_term[τ.<0] .= 0.0
    return  frequancy_term + loss_term # does it have to be multiplied by U ??
end
