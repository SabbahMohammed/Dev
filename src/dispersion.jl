
#=
This file can be used to find betas for air and N2
To find β2 at 800nm for example, use beta(800,2) .If β3 is required use beta(800,3) and so on.
to plot these betas, use plot(eqn, 200:10:1300)
=#



# x is ω
function beta0(x, gas, p, Tk, a)

    # return x/c + 0.5*x*p*Tk_0*(b_1/(1-c_1*x^2) + b_2/(1-c_2*x^2))/(p_0*Tk)/(c) - .5*u^2*c/(a^2*x)
    # return x/c + 0.5*x*(refractive_index(x,gas,p,Tk)-1)/(c) - .5*u^2*c/(a^2*x)
    L = x/c - .5*u^2*c/(a^2*x)
    for i in keys(gas)
        b_1,b_2,c_1,c_2 = dispersion_coeffecient(i)
        L +=  0.5*x*gas[i]*p*Tk_0*(b_1/(1-c_1*x^2) + b_2/(1-c_2*x^2))/(c*p_0*Tk)
    end
    return L
end

function dispersion_coeffecient(gas)
    if gas == :Air
        # constants for Air
        b_1 = 14926.44*10^-8.0
        b_2 = 41807.57*10^-8.0
        c_1 = 19.36*10^-6.0*10^-12.0/(4*pi^2*c^2)
        c_2 = 7.434*10^-3.0*10^-12.0/(4*pi^2*c^2)
    end


    if gas == :N2
        # constants for N2
        b_1 = 39209.95*10^-8.0
        b_2 = 18806.48*10^-8.0
        c_1 = 1146.24*10^-6.0*10^-12.0/(4*pi^2*c^2)
        c_2 = 13.476*10^-3.0*10^-12.0/(4*pi^2*c^2)
    end


    if gas == :Ar
        # constants for Ar
        b_1 = 20332.29*10^-8.0
        b_2 = 34458.31*10^-8.0
        c_1 = 206.12*10^-6.0*10^-12.0/(4*pi^2*c^2)
        c_2 = 8.066*10^-3.0*10^-12.0/(4*pi^2*c^2)
    end

    if gas == :He
        b_1 = 4977.77*10^-8.0
        b_2 = 1856.94*10^-8.0
        c_1 = 28.54*10^-6.0*10^-12.0/(4*pi^2*c^2)
        c_2 = 7.760*10^-3.0*10^-12.0/(4*pi^2*c^2)
    end

    return b_1,b_2,c_1,c_2
end

function refractive_index(gas)
    if gas == :O3
        data = readdlm(raw"C:\Users\ms246\Dropbox (Heriot-Watt University Team)\RES_EPS_Lupo\Projects\Mohammed\phd\ozone absorption spectrum\kk_1nm.txt")
    end
    data = data[1,:]
    freq = readdlm(raw"C:\Users\ms246\Dropbox (Heriot-Watt University Team)\RES_EPS_Lupo\Projects\Mohammed\phd\ozone absorption spectrum\freqdata_1nm.txt")
    freq = freq[:,1]

    spl = Spline1D(freq, data; k=3, bc="nearest", s=0.0)
    return spl
end

function O3_disp(W,p,Tk,a; diff=0)
    # β   = W./c .+ 0.5.*W.*p.*Tk_0.*refractive_index(:O3)./(c.*p_0.*Tk) .- .5*u^2*c/(a^2*W)
    ref = Spline1D(W, refractive_index(:O3); k=3, bc="nearest", s=0.0)
    y(ω) = ω/c + ω*p*Tk_0*ref(ω)/(c*p_0*Tk) +  0.5*ω*p*293*ref(ω)/(c*p_0*Tk) - 0.5*u^2*c/(a^2.0*ω)
    function β(W)
        temp = Array([])
        for i in 1:length(W)
            if W[i] > 2e16
                append!(temp, y(2e16))
            elseif W[i] < 0.5e15
                append!(temp, y(0.5e15))
            else
                append!(temp, y(W[i]))
            end
        end
        return temp
    end
    β0 = Spline1D(W, β(W), ; k=3, bc="nearest", s=0.0)
    β1 = Spline1D(W, derivative(β0, W); k=3, bc="nearest", s=0.0)

    if diff == 0
        return β0(W)
    elseif diff == 1
        return β1(W)
    elseif diff == 2
        β2 = derivative(β1, W)
        return β2
    end
end

function beta1(x, gas, p, Tk, a)
    return diff(beta0(x, gas, p, Tk, a), x)
end

function beta2(x, gas, p, Tk, a)
    return diff(beta1(x, gas, p, Tk, a), x)
end

function wbeta00(ω::Float64, gas, p, Tk, a)
    # if ω > 2e16
    #     y = beta0(x, gas, p, Tk, a)
    #     return Float64(subs(y, x=>2e16))
    # end
    # if ω < 0.5e15
    #     y = beta0(x, gas, p, Tk, a)
    #     return Float64(subs(y, x=>0.5e15))
    # end
    y = beta0(x, gas, p, Tk, a)
    return Float64(subs(y, x=>ω))
end

function wbeta0(W, gas, p, Tk, a)
    y = Vector()
    for i in 1:length(W)
        l = wbeta00(W[i], gas, p, Tk, a)
        append!(y, l)
    end
    return y
end


function lbeta0(λ, gas::Symbol, p, Tk, a)
    ω = (2pi*c)/λ
    y = beta0(x, gas, p, Tk, a)
    return Float64(subs(y, x=>ω))
end

function wbeta1(ω::Float64, gas::Symbol, p, Tk, a)
    if ω > 2e16
        y = beta1(x, gas, p, Tk, a)
        return Float64(subs(y, x=>2e16))
    end
    if ω < 0.5e15
        y = beta1(x, gas, p, Tk, a)
        return Float64(subs(y, x=>0.5e15))
    end
    y = beta1(x, gas, p, Tk, a)
    return Float64(subs(y, x=>ω))
end

function wbeta2(ω::Float64, gas::Symbol, p, Tk, a)
    if ω > 2e16
        y = beta2(x, gas, p, Tk, a)
        return Float64(subs(y, x=>2e16))
    end
    if ω < 0.5e15
        y = beta2(x, gas, p, Tk, a)
        return Float64(subs(y, x=>0.5e15))
    end
    y = beta2(x, gas, p, Tk, a)
    return Float64(subs(y, x=>ω))
end

# function lbeta0(λ, gas::Symbol, p, Tk, a)
#     ω = (2pi*c)/λ
#     y = beta0(x, gas, p, Tk, a)
#     print(y)
#     return Float64(subs(y, x=>ω))
# end

function lbeta2(λ, gas::Symbol, p, Tk, a)
    ω = (2pi*c)/λ
    y = beta2(x, gas, p, Tk, a)
    return Float64(subs(y, (x,ω)))
end


function L_operator(ω::Float64, ω0, gas::Symbol, p, Tk, a)
    s = ω+ω0
    beta00 = beta0(x, gas, p, Tk, a)
    beta11 = beta1(x, gas, p, Tk, a)
    return Float64(subs(beta00, x=>s)) - Float64(subs(beta00, x=>ω0)) - Float64(subs(beta11, x=>ω0))*ω
end

function LE_operator(W, ω0, gas::Dict{Symbol,Float64}, p, Tk, a)
    L = Vector()
    for i in 1:length(W)
        l = LE0_operator(W[i], ω0, gas, p, Tk, a)
        append!(L, l)
    end
    return L
end

function LE0_operator(ω::Float64, ω0, gas::Dict{Symbol,Float64}, p, Tk, a)
    # if ω == 0
    #     return 0.0
    # end

    if ω > 2e16
        beta00 = beta0(x, gas, p, Tk, a)
        beta11 = beta1(x, gas, p, Tk, a)
        return Float64(subs(beta00, x=>2e16))  - Float64(subs(beta11, x=>ω0))*ω
    end
    if ω < 0.5e15
        beta00 = beta0(x, gas, p, Tk, a)
        beta11 = beta1(x, gas, p, Tk, a)
        return Float64(subs(beta00, x=>0.5e15))  - Float64(subs(beta11, x=>ω0))*ω
    end
    beta00 = beta0(x, gas, p, Tk, a)
    beta11 = beta1(x, gas, p, Tk, a)
    return Float64(subs(beta00, x=>ω))  - Float64(subs(beta11, x=>ω0))*ω
end

# function L_operator(ω::Array, ω0, gas::Symbol, p, Tk, a)
#     s = ω.+ω0
#     beta00 = beta0(x, gas, p, Tk, a)
#     beta11 = beta1(x, gas, p, Tk, a)
#     return Float64(subs(beta00, x=>s)) .- Float64(subs(beta00, x=>ω0)) .- Float64(subs(beta11, x=>ω0)).*ω
# end


function beta3(x, gas::Symbol, p, Tk, a)
    return diff(beta2(x, gas, p, Tk, a), x)
end

function lbeta3(λ, gas::Symbol, p, Tk, a)
    ω = (2pi*c)/λ
    y = beta3(x, gas, p, Tk, a)
    return Float64(subs(y, (x,ω)))
end

function βnl(ω, fiber, pulse)
    ω0 = 2pi*c/pulse.λ0
    γ = getGamma(fiber.gas,  fiber.p, fiber.Tk, ω0, fiber.diameter)
    return wbeta0(ω, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - wbeta0(ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) -
           wbeta1(ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2)*(ω - ω0) - ω/ω0*γ*pulse.P0
end

function dis_wave(fiber, pulse; range = (1e15, 1.5e16))
    βnl0(ω) = βnl(ω, fiber, pulse)
    W = find_zero(βnl0, range)
    return 2*pi*c/W
end

function zero_disp(fiber, pulse; range = (1e15, 1.5e16))
    ω0 = 2pi*c/pulse.λ0
    γ = getGamma(fiber.gas,  fiber.p, fiber.Tk, ω0, fiber.diameter)
    var(ω) = wbeta2(ω, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2)#wbeta0(ω, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - wbeta0(ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2) - wbeta1(ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2)*(ω - ω0)
    zero0 = find_zero(var, range)
    return 2*pi*c/zero0
end

function soliton_order(pulse::Pulse, fiber)
    ω0 = 2pi*c/pulse.λ0
    γ = getGamma(fiber.gas,fiber.p,fiber.Tk,ω0,fiber.diameter)
    τ0 = pulse.τ0
    N = sqrt(γ*pulse.P0*τ0^2/abs(wbeta2(ω0,fiber.gas,fiber.p,fiber.Tk,fiber.diameter/2)))
end

function lfiss(pulse, fiber)
    ω0 = 2pi*c/pulse.λ0
    β2 = wbeta2(ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2)
    γ = getGamma(fiber.gas,fiber.p,fiber.Tk,ω0,fiber.diameter)
    return sqrt(pulse.τ0^2/(abs(β2)*γ*pulse.P0))
end

function soliton_period(pulse, fiber)
    ω0 = 2pi*c/pulse.λ0
    β2 = wbeta2(ω0, fiber.gas, fiber.p, fiber.Tk, fiber.diameter/2)
    return pi*pulse.τ0^2/(2*abs(β2))
end
