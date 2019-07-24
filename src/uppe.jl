function uppe(pulse::Pulse, fiber::Fiber, betas::Array{Any,1}; dispersion::Bool = true, thg::Bool = false,  fr::Float64=0.0, shock::Bool=false, plasma::Bool=false, pltype::Symbol=:adk, nsaves::Int64=1000)
    barglyphs="[=> ]"

    ω0 = 2pi*c/pulse.λ0
    # grid
    τ = pulse.τ
    dτ = τ[2] - τ[1]
    convt = length(τ)*dτ
    dW = 2pi/(τ[end]-τ[1])
    W = (0:length(τ)/2)*dW
    convw = 1/dτ
    convt = dτ

    Et = getElectricFieldE(pulse, fiber)# .+ getElectricFieldE(pulse2, fiber)
    Ew0 = rfft(Et)*dτ

    # fftw operators
    planrfft = plan_rfft(Et)
    planirfft = plan_irfft(Ew0,size(Et)[1])
    #
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



    #loss = 0
    #dispersion
    if dispersion == true
        if fiber.gas == :O3
            if length(τ) == 2^14
                return 1im.*O3_func(W,p,Tk,a)
            else
                throw("for O3, grid size should be 2^14")
            end
        end
        L = 1im.*LE_operator(W,ω0, fiber.gas, fiber.p[1], fiber.Tk, fiber.diameter/2) #- fiber.α/2
    else
        L = 0
    end
    #Raman
    if fr != 0.0
        Ht = Htot(τ, fiber.gas,fiber.p[1], fiber.Tk)
        Hw = planrfft*fftshift(Ht)
    else
        Hw = [0.0im]
    end

    #windows
    gw = exp.(-100 .*((W .-  maximum(W)./2.2)./(maximum(W).*0.47)).^40)
    # gw = exp.(-100 .*((W .-  2.8e15)./(maximum(W).*0.11)).^40)


    gt =  exp.(-0.5 .*((τ)./(maximum(τ).*0.8)).^30)

    #some calculation done insde rhs
    wb0 = W.^2.0./(2*c^2.0.*ϵ0.*wbeta0(W,fiber.gas, fiber.p[1], fiber.Tk, fiber.diameter/2))
    wbdiff = wb0[end] - wb0[1]
    wdiff = W[end] - W[1]
    wb = wbdiff/wdiff*W
    χ31 = χ3(fiber.gas, fiber.p[1], fiber.Tk)


    dir = "C:\\Users\\ms246\\OneDrive - Heriot-Watt University\\programming\\julia\\FromGit\\uppe\\"

    open(dir*"L.txt", "w") do io
         writedlm(io, abs2.(L))
    end




    sys = System{Float64}(npts = size(τ)[1], L=L, Et=Et, W=W, τ=τ, dτ=dτ, dW=dW, bar=[1,2,3],
                ω0=ω0, fiber=fiber, pulse=pulse, gw=gw, gt=gt, convw=convw, convt=convt,
                wb=wb, plasma=plasma, thg=thg, planrfft=planrfft, planirfft=planirfft,
                GePl=GePl, Hw=Hw, fr=fr, χ3=χ31, pltype=pltype)



    zs = (0.0, fiber.length)
    Z = LinRange(0, fiber.length, nsaves)
    prob1 = ODEProblem(rhs2, Ew0, zs, (ProgressMeter.Progress(floor(Int, Float64(fiber.length*1000)), dt=2.0, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow),sys))
    # sol = solve(prob1, Tsit5(), reltol=1e-8, abstol=1e-14, force_dtmin=true, saveat=fiber.length/nsaves)
    @time sol = solve(prob1, Tsit5(), reltol=1e-2, abstol=1e-9,
                    saveat=fiber.length/nsaves, alias_u0= true, dense=false)

    res = zeros(Complex{Float64}, size(sol))
    @simd for i in eachindex(sol.t)
        @inbounds res[:,i] = Complex{Float64}.(sol[:,i] .* exp.(-L .* sol.t[i]))#
    end

    rest = zeros(Complex{Float64}, (length(τ), nsaves+1))
    @simd for i in eachindex(Z)
        Ew = res[:,i].*2
        Ew1 = append!(zeros(Complex{Float64}, (length(τ)- length(Ew))), Ew)
        @inbounds rest[:,i] = Complex{Float64}.(ifft(Ew1))
    end

    Ets = zeros(Complex{Float64}, (length(τ), nsaves+1))
    @simd for i in eachindex(Z)
        Ets[:,i] = irfft(res[:,i], size(τ)[1])./sqrt(2/(ϵ0*c*effectiveArea(fiber.diameter)))*1/dτ
    end

    envs = zeros(Float64, (length(τ), nsaves+1))
    @simd for i in eachindex(Z)
        envs[:,i] = abs2.(get_env(Ets[:,i]))
    end

    τ, W, envs, res, Z, L

    # τ, W, res
end



function rhs2(dEw, Ew, para, z)
    sys = para[2]
    # ir = para[3]

    update!(para[1], max(floor(Int, Float64(z*1000)), 1))

    # if z>0.06 && z<0.09
    #     sys.fiber.p = sys.fiber.p+2
    #     sys.L = @. 1im*LE_operator.(W,ω0, fiber.gas, sys.fiber.p, fiber.Tk, fiber.diameter/2, size(τ)[1])
    # end

    if z>0.06
        sys.fiber.gas[:O3] = 1e-3
    end

    sys.Ew_0 .= Ew
    ml!(sys.Ew_0, sys.L, z)

    sys.temp5 .= sys.planirfft*sys.Ew_0
    ift!(sys.Et, sys.temp5 , sys.convw)

    if sys.fr != 0.0
        sys.preRt .=  (sys.planirfft*(sys.Hw.*(sys.planrfft*(sys.Et.^2))))
        sys.Rt .= sys.fr.*sys.Et.*sys.χ3.*ϵ0.*sys.preRt.*sys.dτ
        sys.ret .= (1-sys.fr)*sys.χ3*ϵ0.*sys.Et.^3 .+ sys.Rt
    else
        if sys.thg == true
            sys.ret .= sys.χ3*ϵ0.*sys.Et.^3
        else
            sys.ret .= 0.0
        end
    end

    if sys.plasma == true
        if sys.pltype == :adk
            if sys.fiber.gas == :Air
                GePlN2 = sys.GePl[2]
                GePlO2 = sys.GePl[1]
                geisslerPlasmaUPPE(sys.Et, sys.τ, GePlN2)
                geisslerPlasmaUPPE(sys.Et, sys.τ, GePlO2)
                sys.Pe .= 0.2*GePlO2.Pe + 0.8*GePlN2.Pe
                # geisslerPlasmaUPPE(Et, τ, Pe, Na, Ne, WE, Pe_l, Pe_f, plconstants, temp1, temp2, temp3, temp4)
            else
                geisslerPlasmaUPPE(sys.Et, sys.τ, sys.GePl[1])
                sys.Pe .= sys.GePl[1].Pe
            end
        elseif sys.pltype == :ppt
            println("ppt model is not ready yet")
        end
        sys.ret .= sys.ret .+ sys.Pe
    end
    # else
    #     ml_1!(ret, (Et).^3, χ3(fiber.gas, fiber.p, fiber.Tk)*ϵ0)
    # end

    ml_2!(sys.PNL, sys.wb, sys.planrfft*sys.ret, -1im*sys.convt)

    ml_2!(sys.dEt, sys.planirfft*sys.PNL, sys.gt, sys.convw)

    ml_0!(dEw, sys.planrfft*sys.dEt, sys.L, sys.gw, z, sys.convt)
end






#
# function rhs2(dEw, Ew, para, z)
#     L = para[1]
#     W, τ, ω0, fiber, pulse, gw, gt, convw, convt, Et, ret, PNL, dEt,
#         Pe, Pe0, Na, Ne, WE, Pe_l, Pe_f, wb, plasma, thg, dτ, dW = para[3:end]
#     dτ = τ[2] - τ[1]
#     dW = 2pi/(τ[end]-τ[1])
#     #print(z)
#     update!(para[2], max(floor(Int, Float64(z*1000)), 1))
#
#     #=
#     why the scaling is like this?
#     convw = 1/dτ
#     convt = dτ
#     from the fact that how DFT work (https://en.wikipedia.org/wiki/Discrete_Fourier_transform)
#     it is clear that the algorithm does not include the integral step (ie dt), so we multiply fft
#     by dτ and ifft by df*N which gives 1/dt(the *N is because in the algorithm it divide by N)
#     =#
#
#     Ew = Ew.*exp.(-L.*z)
#
#     Et .= irfft(Ew, length(τ)) .* convw # real.(ifft(cat(conj(Ew[end-1:-1:2]), Ew, dims =1))) .* convw
#
#     if plasma == true
#         Pe .= geisslerPlasmaUPPE(Et, τ, dτ, fiber.Ip, fiber.p, Pe0, Na, Ne, WE, Pe_l, Pe_f)
#         ret .= χ3(fiber.gas, fiber.p, fiber.Tk)*ϵ0.*Et.^3 .+ Pe
#     else
#         ret .= χ3(fiber.gas, fiber.p, fiber.Tk)*ϵ0.*(Et).^3
#     end
#
#     PNL .= -1im.*wb.*rfft(ret) .* convt
#
#     dEt .= irfft(PNL, length(τ)) .* convw .* gt
#
#     # println(z)
#     dEw .= rfft(dEt) .*exp.(L.*z).* convt.* gw
# end
