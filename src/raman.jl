function Nk(J, Bk)
    return exp(-h*c*Bk*J*(J+1)/(KB*Tk))
end

function Hrot(τ, nJ, gas; pre=-1.0536477760095002e12)
    qt = zeros(Float64, size(τ))
    if gas == :N2
        Bk = 1.99*100 #(m^-1) in QSM paper the constant s in cm^-1
        τrot_k = 61e-12
        for j in 1:nJ
            if j % 2 == true
                qk = 1
            else
                qk = 2
            end
            A_k =(Nk(j+2, Bk) - Nk(j, Bk))*qk*(j + 2)*(j + 1)/(2*j + 3)
            qt += @. (exp(-τ/τrot_k)*A_k*sin(4pi*Bk*c*(2*j+3)*τ))
        end
    end
    if gas == :O2
        Bk = 1.44*100 #(m^-1) in QSM paper the constant s in cm^-1
        τrot_k = 68e-12
        for j in 1:nJ
            if j % 2 == true
                qk = 1
            else
                qk = 0
            end
            A_k =(Nk(j+2, Bk) - Nk(j, Bk))*qk*(j + 2)*(j + 1)/(2*j + 3)
            qt += @. (exp(-τ/τrot_k)*A_k*sin(4pi*Bk*c*(2*j+3)*τ))
        end
    end
     # ret[τ<0] = 0
    @. qt[τ<0] = 0
    return qt .* pre
end

function Hvib(τ, nJ, gas; pre=2.8945247923481406e12)
    qt = zeros(Float64, size(τ))

    if gas == :N2
        Bk = 1.99*100
        τvib_k = 10e-12
        Ω_k = 2330.0*100.0 #(m^-1)
        η_k = 0.0173*100.0
        for j in 1:nJ
            if j % 2 == true
                qk = 1
            else
                qk = 2
            end
            M_k = qk*(2*j+1)*exp(-h*c*Bk*j*(j+1)/(KB*Tk))
            ω_k = 2*pi*c*(Ω_k - η_k*j*(j+1))
            @. qt += (exp(-τ/τvib_k)*(M_k*sin(ω_k*τ)))
        end
    end
    if gas == :O2
        Bk = 1.44*100
        τvib_k = 10e-12
        Ω_k = 1556.0*100.0 #(m^-1)
        η_k = 0.0159*100.0
        for j in 1:nJ
            if j % 2 == true
                qk = 1
            else
                qk = 0
            end
            M_k = qk*(2*j+1)*exp(-h*c*Bk*j*(j+1)/(KB*Tk))
            ω_k = 2*pi*c*(Ω_k - η_k*j*(j+1))
            @. qt += (exp(-τ/τvib_k)*(M_k*sin(ω_k*τ)))
        end
    end
    @. qt[τ<0] = 0
    return qt .* pre
end

function Htot(τ, gas, p, Tk; nJ=30)
    Ht = zeros(Float64, size(τ))
    if gas == :Air
        μrot1 = 0.986
        μrot2 = 1
        H1 = μrot1.*Hrot(τ, nJ, :N2) + (1-μrot1).*Hvib(τ, nJ, :O2)
        H2 = μrot2.*Hrot(τ, nJ, :N2) + (1-μrot2).*Hvib(τ, nJ, :O2)
        nrN2 = n_2(:N2,  p, Tk)
        nrO2 = n_2(:O2,  p, Tk)
        nrAir = n_2(:Air,  p, Tk)
        Ht = 0.8*nrN2/nrAir*(μrot1.*Hrot(τ, nJ, :N2) + (1-μrot1).*Hvib(τ, nJ, :N2)) +
             0.2*nrO2/nrAir*(μrot2.*Hrot(τ, nJ, :N2) + (1-μrot2).*Hvib(τ, nJ, :N2))
        return Ht
    elseif gas == :N2
        μrot = 0.986
        μrot.*Hrot(τ, nJ, gas) + (1-μrot).*Hvib(τ, nJ, gas)
    end
end

function getpre(nJ, gas)
    τ = LinRange(-1.0, 1000.0,3*10^6)*1e-12
    Hr = Hrot(τ, nJ, gas, pre=1.0)
    totr = trapz(τ, Hr)
    Hv = Hvib(τ, nJ, gas, pre=1.0)
    totv = trapz(τ, Hv)
    return 1.0/totr, 1.0/totv
end

function sdo(τ; τ1 = 62.5e-15, τ2=77e-15)
    RT = (τ1^2+τ2^2)/(τ1*τ2^2).*exp.(-τ/τ2).*sin.(τ/τ1)
    RT[τ.<0] .= 0.0
    return RT
end

function raman(fr, γ, convt)
    if fr != 0.0
        γ = γ/(1-fr)
        RW = Htot(τ)
        @. qt[τ<0] = 0
        Ht = f2(Ht, τ)
        @. RW[τ<0] = 0
        RW /= trapz(τ, RW)
        RW = shifted_ifft(RW)*convt
    else
        RW = 0.0
    end
    return RW
end

function covRaman(U, RW, convt, convw, fr)
    # println(fr)
    Uw = shifted_ifft(abs2.(U)).*convt
    Uw = shifted_fft(Uw.*RW) .* U .* convw
    Uw = shifted_ifft(Uw) .* convt
    @. Uw =  (1.0 - fr) * Uw + fr* Uw
    return Uw
end
