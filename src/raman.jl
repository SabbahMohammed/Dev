function Nk(J, Bk)
    return exp(-h*c*Bk*J*(J+1)/(KB*Tk))
end

function Hrot(τ, nJ; pre=getpre(τ, nJ)[1])
    ret = fill(0,(size(τ)))
    Bk = 1.99*100 #(m^-1) in QSM paper the constant s in cm^-1
    τrot_k = 61e-12
    for j in 1:nJ
        if j % 2 == true
            qk = 1
        else
            qk = 2
        end
        A_k =(Nk(j+2, Bk) - Nk(j, Bk))*qk*(j + 2)*(j + 1)/(2*j + 3)
        ret += @. (exp(-τ/τrot_k)*A_k*sin(4pi*Bk*c*(2*j+3)*τ))
    end
     # ret[τ<0] = 0
    @. ret[τ<0] = 0
    return ret .* pre
end

function Hvib(τ, nJ; pre=getpre(τ, nJ)[2])
    ret = fill(0,(size(τ)))
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
        @. ret += (exp(-τ/τvib_k)*(M_k*sin(ω_k*τ)))
    end
    @. ret[τ<0] = 0
    return ret .* pre
end

function Htot(τ; μrot = 0.986, nJ=30)
    return μrot.*Hrot(τ, nJ) + (1-μrot).*Hvib(τ, nJ)
end

function getpre(τ, nJ)
    Hr = Hrot(τ, nJ, pre=1.0)
    totr = trapz(τ, Hr)
    Hv = Hvib(τ, nJ, pre=1.0)
    totv = trapz(τ, Hv)
    return 1.0/totr, 1.0/totv
end

function sdo(τ, τ1, τ2)
    RT = (τ1^2+τ2^2)/(τ1*τ2^2).*exp.(-τ/τ2).*sin.(τ/τ1)
    return RT
end

function raman(fr, γ, convt)
    if fr != 0.0
        γ = γ/(1-fr)
        RW = Htot(τ)
        @. ret[τ<0] = 0
        Ht = f2(Ht, τ)
        @. RW[τ<0] = 0
        RW /= integrate(τ, RW)
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
