function gasPlasma(Ip::Float64, npts::Int64)
    ns = sqrt(13.59844*e/Ip)
    wp = Ip/hbar
    cn2 = 2^(2*ns)/(ns*gamma(ns+1)*gamma(ns))
    return GasPlasma(ns=ns,wp=wp,cn2=cn2, npts=npts)
end

function ionRateAdk(Et::Array{Float64,1}, WE::Array{Float64,1}, Ip::Float64, ns::Float64, wp::Float64, cn2::Float64)
    @. WE = e*abs(Et)/sqrt(2*m_e*Ip)
    # for i in eachindex(WE) the same as above
    #     WE[i] = e*abs(Et[i])/sqrt(2*m_e*Ip)
    # end
    @. WE =  wp*cn2*(4wp/WE)^(2ns-1)*exp(-4.0/3.0*wp/WE)
    # for i in eachindex(WE) The same as above
    #     WE[i] = wp*cn2*(4wp/WE[i])^(2ns-1)*exp(-4.0/3.0*wp/WE[i])
    # end
    # @time  WE[notfinite.(WE)] .= 0.0
    for i in eachindex(WE)
        if ~isfinite(WE[i])
            WE[i] = 0.0
        end
    end
end

function adk(Et::Array{Float64,1}, τ::Array{Float64,1}, gas::Gas)
    adk(Et, τ, gas.gp.Na, gas.gp.Ne, gas.gp.WE, gas)
end

function adk(Et::Array{Float64,1}, τ::Array{Float64,1}, Na::Array{Float64,1}, 
            Ne::Array{Float64,1}, WE::Array{Float64,1}, gas::Gas)
    ionRateAdk(Et, gas.gp.WE, gas.Ip, gas.gp.ns, gas.gp.wp, gas.gp.cn2)
    # @time ionRateAdk(Et, gas)
    
    # @time Na .= exp.(-cumul_integrate(τ, WE))
    cumtrapz(τ, WE, gas.gp.temp1)
    
    @. Na = exp(-gas.gp.temp1)
    # @time for i in eachindex(Na) the above is better
    #     Na[i] = exp(-gas.gp.temp1[i])
    # end

    # Ne[notfinite.(Ne) .| (Ne .< 0.0)] .= 0.0
    # Ne .= gas.p*gas.pp*N0.*(1.0.-Na)
    for i in eachindex(Ne)
        Ne[i] = gas.p*gas.pp*N0*(1.0-Na[i])

        if ~isfinite(Ne[i]) | (Ne[i] < 0.0)
            Ne[i] = 0.0
        end
    end
    

    # Na .= gas.p*gas.pp*N0.*Na
    # @time Na[notfinite.(Na) .| (Na .< 0.0)] .= 1.0
    for i in eachindex(Na)
        Na[i] = gas.p*gas.pp*N0*Na[i]

        if ~isfinite(Na[i]) | (Na[i] < 0.0)
            Na[i] = 0.0
        end
    end
end

function geisslerPlasmaUPPE(Et::Array{Float64,1}, τ::Array{Float64,1}, dτ::Float64, gas::Gas)
    adk(Et, τ, gas)
    @. gas.gp.temp2 = gas.gp.Na*gas.gp.WE/Et
    #the above is better. something wrong with passing gas.gp.element to for loop
    # (seems to be a type problem. I need to specify all types to git rid of the problem)
    # @time for i in eachindex(gas.gp.temp1) 
    #     gas.gp.temp2[i] = gas.gp.Na[i]*gas.gp.WE[i]/Et[i]
    # end
    cumtrapz(τ, gas.gp.temp2, gas.gp.Pe_l)
    @. gas.gp.Pe_l .= gas.Ip .* gas.gp.Pe_l


    # gas.gp.Pe_f .= cumul_integrate(τ, gas.gp.Ne.*Et)
    # gas.gp.Pe_f .= e^2/m_e.*cumul_integrate(τ, gas.gp.Pe_f)
    @. gas.gp.temp3 = gas.gp.Ne.*Et
    cumtrapz(τ, gas.gp.temp3, gas.gp.temp4)
    cumtrapz(τ, gas.gp.temp4, gas.gp.Pe_f)


    @. gas.gp.Pe = e^2/m_e*gas.gp.Pe_f + gas.gp.Pe_l
    @. gas.gp.Pe[notfinite.(gas.gp.Pe)] = 0.0
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


function ppt(Et, ir::IonizationRatePPT; tol::Float64 = 1e-3, rcycle::Bool = true)
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


for gnlse
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
    



   # old plasma gas mixture for air
    # if fiber.gas == :Air
    #
    #     O2_Ip = 12.0697*e
    #     N2_Ip = 15.581*e
    #
    #     nsO2 = sqrt(13.59844*e/O2_Ip)
    #     wpO2 = O2_Ip/hbar
    #     cn2O2 = 2^(2*nsO2)/(nsO2*gamma(nsO2+1)*gamma(nsO2))
    #
    #     nsN2 = sqrt(13.59844*e/N2_Ip)
    #     wpN2 = N2_Ip/hbar
    #     cn2N2 = 2^(2*nsN2)/(nsN2*gamma(nsN2+1)*gamma(nsN2))
    #
    #     plconstantsO2 = [O2_Ip, nsO2, wpO2, cn2O2, fiber.p]
    #     plconstantsN2 = [N2_Ip, nsN2, wpN2, cn2N2, fiber.p[1]]
    #
    #     GePlO2 = geisslerPlasmaUPPE{Float64}(npts = size(τ)[1], p = fiber.p[1], plconstants = plconstantsO2)
    #     GePlN2 = geisslerPlasmaUPPE{Float64}(npts = size(τ)[1], p = fiber.p[1], plconstants = plconstantsN2)
    #     GePl = [GePlO2, GePlN2]
    # else
    #     Ip = fiber.Ip*e
    #     ns = sqrt(13.59844*e/Ip)
    #     wp = Ip/hbar
    #     cn2 = 2^(2*ns)/(ns*gamma(ns+1)*gamma(ns))
    #     plconstants = [Ip, ns, wp, cn2, fiber.p[1]]
    #     ir = IonizationRatePPT(fiber.Ip, 3, ω0; Z=1)
    #     GePl = [geisslerPlasmaUPPE{Float64}(npts = size(τ)[1],plconstants = plconstants, pltype=pltype, p=fiber.p[1])]
    # end
    # if plasma == true
    #     Ip = fiber.Ip*e
    #     ns = sqrt(13.59844*e/Ip)
    #     wp = Ip/hbar
    #     cn2 = 2^(2*ns)/(ns*gamma(ns+1)*gamma(ns))
    #     plconstants = [Ip, ns, wp, cn2, fiber.p[1]]
    #     ir = IonizationRatePPT(fiber.Ip, 3, g.ω0; Z=1)
    #     GePl = [geisslerPlasmaUPPE{Float64}(npts = size(τ)[1],plconstants = plconstants, pltype=pltype, p=fiber.p[1])]
    # else
    #     GePl =   [geisslerPlasmaUPPE{Float64}( npts=g.npts, plconstants=[],  p=gases[1].p)]
    # end

    
# function geisslerPlasmaUPPE(Et, τ, gas, pltype)
#     gePl = gas.gasplasma

#     Na, Ne, WE = adk(Et, τ, gas.Ip, gas.p*gas.pp)
#     gePl.WE .= WE
#     gePl.Na .= Na
#     gePl.Ne .= Ne
#     gePl.Pe_l .= gas.Ip .* cumul_integrate(τ, gePl.Na.*gePl.WE./Et)
#     gePl.Pe_f .= cumul_integrate(τ, gePl.Ne.*Et)
#     gePl.Pe_f .= e^2/m_e.*cumul_integrate(τ, gePl.Pe_f)
#     gePl.Pe .= gePl.Pe_f .+ gePl.Pe_l
#     gePl.Pe[notfinite.(gePl.Pe)] .= 0.0

#     # @. gePl.temp3 = gePl.Na*gePl.WE/Et
    

#     # # @time for i  in eachindex(temp3)
#     # #     temp3[i] = Na[i]*WE[i]/Et[i]
#     # # end
#     # cumtrapz(τ, gePl.temp3, gePl.Pe_l)
#     # println(gePl.Pe_l)
    
#     # @. gePl.Pe_l = gas.Ip*gePl.Pe_l

    
#     # # @time for i in eachindex(Pe_l)
#     # #     Pe_l[i] = Ip*Pe_l[i]
#     # # end
#     # @. gePl.temp4 = gePl.Ne*Et
#     # # @time for i in eachindex(temp4)
#     # #     temp4[i] = Ne[i]*Et[i]
#     # # end
#     # cumtrapz(τ, gePl.temp4, gePl.temp1)
#     # cumtrapz(τ, gePl.temp1, gePl.Pe_f)
#     # for i in eachindex(gePl.Pe)
#     #     gePl.Pe[i] = e^2/m_e*gePl.Pe_f[i] + gePl.Pe_l[i]
#     #     if ~isfinite(gePl.Pe[i])
#     #         gePl.Pe[i] = 0.0
#     #     end
#     # end
# end



#=
all below is for uppe
=#
# function ionRateAdk(Et, Ip, gePl)
#     @. gePl.WE = e*abs(Et)/sqrt(2*m_e*Ip)
#     @. gePl.WE =  gePl.wp*gePl.cn2*(4*gePl.wp/gePl.WE)^(2*gePl.ns-1)*exp(-4.0/3.0*gePl.wp/gePl.WE)
#     # @time  WE[notfinite.(WE)] .= 0.0
#     for i in eachindex(gePl.WE)
#         if ~isfinite(gePl.WE[i])
#             gePl.WE[i] = 0.0
#         end
#     end
# end

# function ionRateAdk(Et, WE, plconstants)
#     Ip, ns, wp, cn2, p = plconstants
#     @. WE = e*abs(Et)/sqrt(2*m_e*Ip)
#     @. WE =  wp*cn2*(4wp/WE)^(2ns-1)*exp(-4.0/3.0*wp/WE)
#     # @time  WE[notfinite.(WE)] .= 0.0
#     for i in eachindex(WE)
#         if ~isfinite(WE[i])
#             WE[i] = 0.0
#         end
#     end
# end

# initial gas mixture try (did not work)
# function ionRate(Et, τ, gePl, gas,pltype)
#     if pltype == :adk
#         ionRateAdk(Et, gas.Ip, gePl)
#     elseif pltype == :ppt
#         ppt.(Et, gePl.ir::IonizationRatePPT)
#     end
#     # p = gePl.plconstants[5]
#     cumtrapz(τ, gePl.WE, gePl.temp2)
#     # Na .= exp.(-temp2)
#     for i in eachindex(gePl.Na)
#         gePl.Na[i] = exp(-gePl.temp2[i])
#     end
#     for i in eachindex(gePl.Ne)
#         # gePl.Ne[i] = p*N0*(1.0-gePl.Na[i]
#         gePl.Ne[i] = gas.p*gas.pp*N0*(1.0-gePl.Na[i])

#         if ~isfinite(gePl.Ne[i]) | (gePl.Ne[i] < 0.0)
#             gePl.Ne[i] = 0.0
#         end
#     end
#     for i in eachindex(gePl.Na)
#         gePl.Na[i] = gas.p*gas.pp*N0*gePl.Na[i]
#         if ~isfinite(gePl.Na[i]) | (gePl.Na[i] < 0.0)
#             gePl.Na[i] = 0.0
#         end
#     end
# end

# function geisslerPlasmaUPPE(Et, τ, dτ, Ip, p, Pe0, Na, Ne, WE, Pe_l, Pe_f)
#     Na, Ne, WE = adk(Et, τ, Ip, p)
#     Pe_l .= Ip .* cumul_integrate(τ, Na.*WE./Et)
#     Pe_f .= cumul_integrate(τ, Ne.*Et)
#     Pe_f .= e^2/m_e.*cumul_integrate(τ, Pe_f)
#     Pe0 .= Pe_f .+ Pe_l
#     Pe0[notfinite.(Pe0)] .= 0.0
#     return Pe0
# end

# function ionRateAdk(Et, Ip)

#     ns = sqrt(13.59844*e/Ip)
#     wp = Ip/hbar
#     cn2 = 2^(2*ns)/(ns*gamma(ns+1)*gamma(ns))
#     wt = @. e*abs(Et)/sqrt(2*m_e*Ip)
#     WE = @. wp*cn2*(4wp/wt)^(2ns-1)*exp(-4.0/3.0*wp/wt)
#     WE[notfinite.(WE)] .= 0.0
#     WE
# end

# function adk(Et, τ, Ip, p)
#     WE = ionRateAdk(Et, Ip)
#     Na = exp.(-cumul_integrate(τ, WE))
#     Ne = p*N0.*(1.0.-Na)
#     Na = p*N0.*Na
#     Na[notfinite.(Na) .| (Na .< 0.0)] .= 1.0
#     Ne[notfinite.(Ne) .| (Ne .< 0.0)] .= 0.0
#     [Na, Ne, WE]
# end

# function ionRateAdk(Et, gas)
#     @. gas.gp.wt =  e*abs(Et)/sqrt(2*m_e*gas.Ip)
#     @. gas.gp.WE =  gas.gp.wp*gas.gp.cn2*(4*gas.gp.wp/gas.gp.wt)^(2*gas.gp.ns-1)*exp(-4.0/3.0*gas.gp.wp/gas.gp.wt)
#     gas.gp.WE[notfinite.(gas.gp.WE)] .= 0.0
# end
