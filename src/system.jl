
@with_kw struct System{T}
    foo::Int = 6
    npts::Float64
    L::Array{Complex{T},1}
    Et::Array{T,1}
    W::Array{T,1}
    τ::Array{T,1}
    dτ::Float64
    dW::Float64
    ω0::Float64
    fiber::Fiber{T, Symbol}
    pulse::Pulse{T}
    χ3::Float64
    gw::Array{T,1}
    gt::Array{T,1}
    convw::Float64
    convt::Float64
    wb::Array{T,1}
    plasma::Bool
    thg::Bool
    planrfft::FFTW.rFFTWPlan{Float64,-1,false,1}
    planirfft::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,1},Float64}
    GePl::Array{geisslerPlasmaUPPE{T},1}
    Hw::Array{Complex{Float64},1}
    fr::Float64
    ret::Array{T,1} = zeros(Float64, npts)#
    Pe0::Array{T,1} = zeros(Float64, npts)
    Pe::Array{T,1}= zeros(Float64, npts)
    Na::Array{T,1}= zeros(Float64, npts)
    Ne::Array{T,1}= zeros(Float64, npts)
    WE::Array{T,1}= zeros(Float64, npts)
    temp1::Array{T,1} = zeros(Float64, npts)
    temp2::Array{T,1} = zeros(Float64, npts)
    temp3::Array{T,1} = zeros(Float64, npts)
    temp4::Array{T,1} = zeros(Float64, npts)
    temp5::Array{T,1} = zeros(Float64, npts)
    Pe_l::Array{T,1} = zeros(Float64, npts)
    Pe_f::Array{T,1} = zeros(Float64, npts)
    PNL::Array{Complex{T},1} = zeros(Complex{Float64}, Int(npts/2+1))
    dEt::Array{T,1} = zeros(Float64, npts)
    Ew_0::Array{Complex{T},1} = zeros(Complex{Float64}, Int(npts/2+1))
    Ew2::Array{Complex{T},1} = zeros(Complex{Float64}, Int(npts/2+1))
    Rt::Array{T,1} = zeros(Float64, npts)
    preRt::Array{T,1} = zeros(Float64, npts)
    pltype::Symbol
end
