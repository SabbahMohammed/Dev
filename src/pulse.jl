struct Pulse{T}
    data #pulse envelope
    FWHM::T
    λ0::T
    τ::Array{T,1}
    P0::T # peak power
    τ0::T
    energy::T
end

function createTimeGrid(time_window; npts = 2^14)
    dt = time_window/npts
    τ = LinRange(-time_window/2, time_window/2, npts)
    τ = collect(τ)
end

function gaussianPulse(;τ= τ, energy=energy, FWHM=FWHM, λ0=λ0)
    τ0 = FWHM/1.664
    P0 = 0.94*energy/FWHM
    res0 = sqrt(P0).*sqrt.((exp.(-(τ.^2.0) ./ (τ0^2.0))))
    res::Pulse = Pulse{Float64}(res0, FWHM, λ0, τ, P0, τ0, energy)
    res
end

function chirpedGaussian(;τ= τ, energy=energy, FWHM=FWHM, λ0=λ0, C=C)
    τ0 = FWHM/1.664
    P0 = 0.94*energy/FWHM
    res0 = sqrt(P0).*sqrt.((exp.(-(1+1im*C)*(τ.^2.0) ./ (τ0^2.0))))
    res::Pulse = Pulse{Float64}(res0, FWHM, λ0, τ, P0, energy)
    res
end

function sechPulse(;τ= τ, energy=energy, FWHM=FWHM, λ0=λ0)
    τ0 = FWHM/1.76
    P0 = energy/(2*τ0)
    res0 = sqrt(P0).*sech.(τ/τ0)
    res::Pulse = Pulse{Float64}(res0, FWHM, λ0, τ, P0, τ0, energy)
    res
end

function pulseFromFile(file::String; timeFactor=1e-15, npts=2^14, energy=1e-6, phase_sign = 1.0)
    data = readdlm(file, Float64)
    # if length(data[1,:]) < 20
    #     time = data[1,:]
    #     intensity = data[2,:]
    #     phase = data[3,:]
    # else
    #     time = data[:1,]
    #     intensity = data[:2,]
    #     phase = data[:3,]
    # end
    timep = data[:,1]
    intensity = data[:,2]
    phase = data[:,3]
    timep *= timeFactor
    timeGrid = LinRange(-1e-12,1e-12,npts)
    # print(time)
    spl_intensity = Spline1D(timep, intensity)
    spl_phase = Spline1D(timep, phase)
    generated_intensity = spl_intensity(timeGrid)
    generated_phase = spl_phase(timeGrid)
    generated_intensity[generated_intensity.<0] .= 0
    xprofile = sqrt.(generated_intensity) .* exp.(phase_sign*1im .* generated_phase)
    norm = trapz(timeGrid, abs2.(xprofile))
    res = xprofile.*sqrt(energy)./sqrt(norm)
    P0 = maximum(abs2.(res))
    fwhm = FWHM(timeGrid, abs2.(res))
    τ0 = fwhm/1.664
    pulse = Pulse{Float64}(res, fwhm, 800e-9, timeGrid, P0, τ0, energy)
    return pulse
end

function FWHM(pulse::Pulse)
    # for this function to be accurte, the number of points in the time grid must be >= 2^14
    a = (pulse.data).^2
    hm = maximum(a)/2.0
    out = pulse.τ[findlast(x -> x > hm, a)] - pulse.τ[findfirst(x -> x > hm, a)]
    out
end

function FWHM(τ, data)
    # for this function to be accurte, the number of points in the time grid must be >= 2^14
    a = abs2.(data)
    hm = maximum(a)/2.0
    out = τ[findlast(x -> x > hm, a)] - τ[findfirst(x -> x > hm, a)]
    out
end

function Waist(τ, data)
    # for this function to be accurte, the number of points in the time grid must be >= 2^14
    a = abs2.(data)
    hm = maximum(a)/2.718281828459045
    out = τ[findlast(x -> x > hm, a)] - τ[findfirst(x -> x > hm, a)]
    out
end

function getElectricField(pulse::Pulse, fiber) #avergae electric field
    ω0 = 2pi*c/pulse.λ0
    Et = @. pulse.data * exp(-1im * ω0 * pulse.τ)
    Et = real(Et)
    Et = sqrt(2/(ϵ0*c*effectiveArea(fiber.diameter)))*Et
    return Et
end
