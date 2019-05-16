#for gnlse
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

#for gnlse
function adk(Et, τ, Ip, p)
    WE = ionRateAdk(Et, Ip)
    Na = exp.(-cumul_integrate(τ, WE))
    Ne = p*N0.*(1.0.-Na)
    Na = p*N0.*Na
    Na[notfinite.(Na) .| (Na .< 0.0)] .= 1.0
    # for i in eachindex(Na)
    Ne[notfinite.(Ne) .| (Ne .< 0.0)] .= 0.0
    [Na, Ne, WE]
end

#for gnlse
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

function ionRateAdk(Et, WE, plconstants)
    Ip, ns, wp, cn2, p = plconstants
    @. WE = e*abs(Et)/sqrt(2*m_e*Ip)
    @. WE =  wp*cn2*(4wp/WE)^(2ns-1)*exp(-4.0/3.0*wp/WE)
    # @time  WE[notfinite.(WE)] .= 0.0
    for i in eachindex(WE)
        if ~isfinite(WE[i])
            WE[i] = 0.0
        end
    end
end

function ionRate(Et, τ, GePl)
    if GePl.pltype == :adk
        ionRateAdk(Et, GePl.WE, GePl.plconstants)
    elseif GePl.pltype == :ppt
        ppt(Et, τ, ir::IonizationRatePPT)
    end
    p = GePl.plconstants[5]
    cumtrapz(τ, GePl.WE, GePl.temp2)
    # Na .= exp.(-temp2)
    for i in eachindex(GePl.Na)
        GePl.Na[i] = exp(-GePl.temp2[i])
    end
    for i in eachindex(GePl.Ne)
        GePl.Ne[i] = p*N0*(1.0-GePl.Na[i])
        if ~isfinite(GePl.Ne[i]) | (GePl.Ne[i] < 0.0)
            GePl.Ne[i] = 0.0
        end
    end
    for i in eachindex(GePl.Na)
        GePl.Na[i] = p*N0*GePl.Na[i]
        if ~isfinite(GePl.Na[i]) | (GePl.Na[i] < 0.0)
            GePl.Na[i] = 0.0
        end
    end
end

@with_kw struct geisslerPlasmaUPPE{T}
    foo::Int = 6
    npts::Float64
    plconstants
    Pe::Array{T,1} = zeros(Float64, npts)
    Na::Array{T,1} = zeros(Float64, npts)
    Ne::Array{T,1} = zeros(Float64, npts)
    WE::Array{T,1} = zeros(Float64, npts)
    Pe_l::Array{T,1} = zeros(Float64, npts)
    Pe_f::Array{T,1} = zeros(Float64, npts)
    temp1::Array{T,1} = zeros(Float64, npts)
    temp2::Array{T,1} = zeros(Float64, npts)
    temp3::Array{T,1} = zeros(Float64, npts)
    temp4::Array{T,1} = zeros(Float64, npts)
    pltype::Symbol = :adk
end




function geisslerPlasmaUPPE(Et, τ, GePl)
    ionRate(Et, τ, GePl)

    Ip = GePl.plconstants[1]
    @. GePl.temp3 = GePl.Na*GePl.WE/Et
    # @time for i  in eachindex(temp3)
    #     temp3[i] = Na[i]*WE[i]/Et[i]
    # end
    cumtrapz(τ, GePl.temp3, GePl.Pe_l)
    @. GePl.Pe_l = Ip*GePl.Pe_l


    # @time for i in eachindex(Pe_l)
    #     Pe_l[i] = Ip*Pe_l[i]
    # end
    @. GePl.temp4 = GePl.Ne*Et
    # @time for i in eachindex(temp4)
    #     temp4[i] = Ne[i]*Et[i]
    # end
    cumtrapz(τ, GePl.temp4, GePl.temp1)
    cumtrapz(τ, GePl.temp1, GePl.Pe_f)
    for i in eachindex(GePl.Pe)
        GePl.Pe[i] = e^2/m_e*GePl.Pe_f[i] + GePl.Pe_l[i]
        if ~isfinite(GePl.Pe[i])
            GePl.Pe[i] = 0.0
        end
    end
end

struct IonizationRatePPT
    ip::Float64           # ionization potential [J]
    Ui_a::Float64         # ionization potential [au]
    E0_a::Float64         # atomic field strength [au]
    ns::Float64           # effective principle quantum number
    ls::Float64           # effective azimuthal quantum number
    Cnl2::Float64         #
    w_a::Float64          # frequency [au]
    l::Int                # azimuthal quantum number
end


function IonizationRatePPT(Ip, n, w; Z=1)
    l = n-1
    ip = Ip*e
    Ui_a = ip/eh
    E0_a = (2Ui_a)^(3/2)
    ns = Z/sqrt(2Ui_a)
    ls = ns - 1
    Cnl2 = 2^(2ns)/(ns*gamma(ns + ls + 1)*gamma(ns - ls))
    w_a = w*a_t
    IonizationRatePPT(ip, Ui_a, E0_a, ns, ls, Cnl2, w_a, l)
end


function ppt(Et, τ, ir::IonizationRatePPT; tol::Float64 = 1e-3, rcycle::Bool = true)
    E_a = abs(E)/a_E
    Eratio = E_a/ir.E0_a
    #if Eratio > 0.5
    #    error("electric field intensity beyond PPT validity: $Eratio")
    #end
    g = ir.w_a*sqrt(2ir.Ui_a)/E_a
    g2 = g*g
    beta = 2g/sqrt(1 + g2)
    alpha = 2(asinh(g) - g/sqrt(1 + g2))
    Up_a = E_a*E_a/(4ir.w_a*ir.w_a)
    Uit_a = ir.Ui_a + Up_a
    v = Uit_a/ir.w_a
    G = 3/(2g)*((1 + 1/(2g2))*asinh(g) - sqrt(1 + g2)/(2g))
    ret = 0.0
    div = 0.0
    for m=-ir.l:ir.l
        div += 1
        am = abs(m)
        flm = ((2ir.l + 1)*factorial(ir.l + am)
               / (2^am*factorial(am)*factorial(ir.l - am)))
        Am = 4/(sqrt(3pi)*factorial(am))*g2/(1 + g2)
        lret = sqrt(3/(2pi))*ir.Cnl2*flm*ir.Ui_a
        lret *= (2ir.E0_a/(E_a*sqrt(1 + g2)))^(2ir.ns - am - 3/2)
        lret *= Am*exp(-2ir.E0_a*G/(3E_a))
        # remove cycle average # TODO: check this
        if rcycle
            lret *= sqrt(pi*ir.E0_a/(3E_a))
        end
        k = ceil(v)
        sum = 0.0
        dsum = 0.0
        while true
            diff = k - v
            dsum = exp(-alpha*diff)*phi(m, sqrt(beta*diff))
            sum += dsum
            k += 1.0
            dsum/sum < tol && break
        end
        lret *= sum
        ret += lret
    end
    ret/(a_t*div)
end
