@with_kw struct Fiber{T, S}
   foo::Int = 6
   diameter::T#fiber diameter
   length::T# fiber length
   gas::S# gas used
   p::T# pressure
   Tk::T# temperture
   Ip::T# ionization energy divided by e. For N2 for example Ip=15.58
   α::T# attenuation Need modification to accept array
end

function systemInfo(pulse::Pulse, fiber::Fiber)
    dict = Dict("FWHM" => pulse.FWHM,
    "wavelength" => pulse.λ0,
    "P0" => pulse.P0,
    "natural_pulse_width" => pulse.τ0,
    "energy" => pulse.energy,
    "gas" => fiber.gas,
    "pressure" => fiber.p,
    "radius" => fiber.diameter/2)
    return JSON.json(dict)

end

function effectiveArea(diameter)
    return 1.5*(diameter/2)^2
end

function rNGas(p, Tk)
    T0 = 273.15
    p0 = 1.01325
    return p*T0/(p0*Tk)
end

function χ3(gas::Symbol, p, Tk)
    #=
    Chi3 of a 'gas' [m^2/V^2], mostly from Lehmeier 1985
    pressure: P [bar]
    temperature: T [K]
    =#
    Chi3He = 3.43e-28
    Chi3N2 = 21.1*Chi3He
    Chi3Ne = 1.8*Chi3He
    Chi3Ar = 23.5*Chi3He
    Chi3Kr = 64.0*Chi3He
    Chi3Xe = 188.2*Chi3He

    if gas == :N2
        return rNGas(p, Tk)*4*Chi3N2
    end
    if gas == :He
        return rNGas(p, Tk)*4*Chi3He
    end
    if gas == :Ne
        return rNGas(p, Tk)*4*Chi3Ne
    end
    if gas == :Ar
        return rNGas(p, Tk)*4*Chi3Ar
    end
    if gas == :Kr
        return rNGas(p, Tk)*4*Chi3Kr
    end
    if gas == :Xe
        return rNGas(p, Tk)*4*Chi3Xe
    end
end


function n_2(gas::Symbol,  p, Tk)
    if gas == :air
        return 1e-23
    end
    res = 3.0*χ3(gas, p, Tk)/(4.0*ϵ0*c)
    return res
end

function getGamma(gas::Symbol,  p, Tk, ω0, diameter)
    return n_2(gas, p, Tk)*ω0/(c*effectiveArea(diameter))
end
