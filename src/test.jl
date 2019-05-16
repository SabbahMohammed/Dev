function get_env(ef)
    ew = fft(ef)
    ew[Int(size(ew)[1]/2+1):end] .= 0.0
    ep = 2.0.*ifft(ew)
    ep
end

Ets = zeros(Complex{Float64}, (length(τ), 100+1))
@simd for i in eachindex(Z)
    Ets[:,i] = irfft(sol[:,i], size(τ)[1])./sqrt(2/(ϵ0*c*effectiveArea(fiber.diameter)))*1/dτ
end
Ets = irfft(sol[:,1], size(τ)(1))
