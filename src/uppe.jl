function uppe(pulse::Pulse, fiber::Fiber, grid::Grid, gases::Array{Gas{Float64, Symbol}}, betas::Array{Any,1}; dispersion::Bool = true, thg::Bool = false,  fr::Float64=0.0, shock::Bool=false, plasma::Bool=false, pltype::Symbol=:adk, nsaves::Int64=1000)
    barglyphs="[=> ]"

    g = grid

    totalp = 0
    for i in 1:length(gases)
        totalp += gases[i].pp
    end

    if totalp != 1
        throw("The total pressure != 1!")
    end

    Et = getElectricFieldE(pulse, fiber)# .+ getElectricFieldE(pulse2, fiber)
    Ew0 = rfft(Et)*g.dτ

    # fftw operators
    planrfft = plan_rfft(Et)
    planirfft = plan_irfft(Ew0,g.npts)
 

    #loss = 0
    #dispersion
    if dispersion == true
        spl = splbeta(g.W, gases, fiber.diameter/2; diff=0)
        L = 1im.*LE_operator(g.W,g.ω0, gases, fiber.diameter/2, spl) #- fiber.α/2
        # data = readdlm("L.txt")
        # L = tryparse.(Complex{Float64}, data[:,2] .* data[:,3])
    else
        L = 0
    end
    #Raman
    if fr != 0.0
        Ht = Htot(g.τ, fiber.gas,fiber.p[1], fiber.Tk)
        Hw = planrfft*fftshift(Ht)
    else
        Hw = [0.0im]
    end

    #windows
    gw = exp.(-100 .*((g.W .-  maximum(g.W)./2.2)./(maximum(g.W).*0.47)).^40)
    gt =  exp.(-0.5 .*((g.τ)./(maximum(g.τ).*0.8)).^30)

    #some calculation done insde rhs
    wb0 = g.W.^2.0./(2*c^2.0.*ϵ0.*beta(g.W, gases, fiber.diameter/2, spl))
    wbdiff = wb0[end] - wb0[1]
    wdiff = g.W[end] - g.W[1]
    wb = wbdiff/wdiff*g.W
    # wb = readdlm("wb.txt")[:,1]
    χ3 = 0
    for i in 1:length(gases)
        χ3 += gases[i].χ3
    end

    sys = System{Float64}(L=L, Et=Et, g=g,
                fiber=fiber, pulse=pulse, gw=gw, gt=gt,wb=wb, plasma=plasma,
                thg=thg, planrfft=planrfft, planirfft=planirfft,
                gases=gases, Hw=Hw, fr=fr, χ3=χ3, pltype=pltype)

    zs = (0.0, fiber.length)
    Z = LinRange(0, fiber.length, nsaves)
    prob1 = ODEProblem(rhs2, Ew0, zs, (ProgressMeter.Progress(floor(Int, Float64(fiber.length*1000)),
                            dt=2.0, barglyphs=BarGlyphs("[=> ]"), barlen=50, color=:yellow),sys))

    @time sol = solve(prob1, Tsit5(), reltol=1e-3, abstol=1e-9,
                    saveat=fiber.length/nsaves, alias_u0= true, dense=false)

    res = zeros(Complex{Float64}, size(sol))
    @simd for i in eachindex(sol.t)
        @inbounds res[:,i] = Complex{Float64}.(sol[:,i] .* exp.(-L .* sol.t[i]))#
    end

    rest = zeros(Complex{Float64}, (g.npts, nsaves+1))
    @simd for i in eachindex(Z)
        Ew = res[:,i].*2
        Ew1 = append!(zeros(Complex{Float64}, (g.npts- g.wnpts)), Ew)
        @inbounds rest[:,i] = Complex{Float64}.(ifft(Ew1))
    end

    Ets = zeros(Complex{Float64}, (g.npts, nsaves+1))
    @simd for i in eachindex(Z)
        Ets[:,i] = irfft(res[:,i], g.npts)./sqrt(2/(ϵ0*c*effectiveArea(fiber.diameter)))*g.convw
    end

    envs = zeros(Float64, (g.npts, nsaves+1))
    @simd for i in eachindex(Z)
        envs[:,i] = abs2.(get_env(Ets[:,i]))
    end

    g.τ, g.W, envs, res, Z, L, sys

    # τ, W, res
end



function rhs2(dEw, Ew, para, z)
    sys = para[2]
    # ir = para[3]

    update!(para[1], max(floor(Int, Float64(z*1000)), 1))

    sys.Ew_0 .= Ew
    ml!(sys.Ew_0, sys.L, z)
    sys.temp5 .= sys.planirfft*sys.Ew_0
    ift!(sys.Et, sys.temp5 , sys.g.convw)

    if sys.fr != 0.0
        sys.preRt .=  (sys.planirfft*(sys.Hw.*(sys.planrfft*(sys.Et.^2))))
        sys.Rt .= sys.fr.*sys.Et.*sys.χ3.*ϵ0.*sys.preRt.*sys.g.dτ
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
            sys.Pe .= zeros(Float64, sys.g.npts)
            for i in eachindex(sys.gases)
                geisslerPlasmaUPPE(sys.Et, sys.g.τ, sys.g.dτ, sys.gases[i])
                sys.Pe .+= sys.gases[i].gp.Pe              
            end
        elseif sys.pltype == :ppt
            println("ppt model is not ready yet")
        end
        sys.ret .= sys.ret .+ sys.Pe
    end

    ml_2!(sys.PNL, sys.wb, sys.planrfft*sys.ret, -1im*sys.g.convt)

    ml_2!(sys.dEt, sys.planirfft*sys.PNL, sys.gt, sys.g.convw)

    ml_0!(dEw, sys.planrfft*sys.dEt, sys.L, sys.gw, z, sys.g.convt)
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
