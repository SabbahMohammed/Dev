@with_kw struct Fiber{T, S}
   foo::Int = 6
   diameter::T#fiber diameter
   length::T# fiber length
   gas::Dict{S,T}# gas used
   p::Array{T,1}# pressure
   Tk::T# temperture
   Ip::T# ionization energy divided by e. For N2 for example Ip=15.58
   α# attenuation Need modification to accept array
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

function χ30(gas, p, Tk)
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
    Chi3O2 = 20.93*Chi3He #estimated from n2 from https://www.osapublishing.org/view_article.cfm?gotourl=https%3A%2F%2Fwww%2Eosapublishing%2Eorg%2FDirectPDFAccess%2FE77D9C07-DA73-08A6-1DBCC7FF6E9E255C_34926%2Fjosab-14-3-650%2Epdf%3Fda%3D1%26id%3D34926%26seq%3D0%26mobile%3Dno&org=Heriot-Watt%20University%20Library
    Chi3Air = Chi3O2*0.20946 + Chi3N2*(1-0.20946-0.009340) + Chi3Ar*0.009340 #estimated from N2 and O2

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
    if gas == :O2
        return rNGas(p, Tk)*4*Chi3O2
    end
    if gas == :Air
        return rNGas(p, Tk)*4*Chi3Air
    end
end

function χ3(gas, p, Tk)
    sus = 0
    for i in keys(gas)
        sus += χ30(i, gas[i]*p, Tk)
    end
    return sus
end



function n_2(gas::Symbol,  p, Tk)
    res = 3.0*χ3(gas, p, Tk)/(4.0*ϵ0*c)
    return res
end

function getGamma(gas::Symbol,  p, Tk, ω0, diameter)
    return n_2(gas, p, Tk)*ω0/(c*effectiveArea(diameter))
end
