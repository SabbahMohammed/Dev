
@with_kw struct System{T}
    foo::Int = 6
    L::Array{Complex{T},1}
    Et::Array{T,1}
    fiber::Fiber{T, Symbol}
    pulse::Pulse{T}
    g::Grid{T}
    gases::Array{Gas,1}
    χ3::Float64
    gw::Array{T,1}
    gt::Array{T,1}
    wb::Array{T,1}
    plasma::Bool
    thg::Bool
    planrfft::FFTW.rFFTWPlan{Float64,-1,false,1}
    planirfft::AbstractFFTs.ScaledPlan{Complex{Float64},FFTW.rFFTWPlan{Complex{Float64},1,false,1},Float64}
    Hw::Array{Complex{Float64},1}
    fr::Float64
    ret::Array{T,1} = zeros(Float64, g.npts)#
    Pe0::Array{T,1} = zeros(Float64, g.npts)
    Pe::Array{T,1}= zeros(Float64, g.npts)
    PeTemp::Array{T,1}= zeros(Float64, g.npts)
    Na::Array{T,1}= zeros(Float64, g.npts)
    Ne::Array{T,1}= zeros(Float64, g.npts)
    WE::Array{T,1}= zeros(Float64, g.npts)
    temp1::Array{T,1} = zeros(Float64, g.npts)
    temp2::Array{T,1} = zeros(Float64, g.npts)
    temp3::Array{T,1} = zeros(Float64, g.npts)
    temp4::Array{T,1} = zeros(Float64, g.npts)
    temp5::Array{T,1} = zeros(Float64, g.npts)
    Pe_l::Array{T,1} = zeros(Float64, g.npts)
    Pe_f::Array{T,1} = zeros(Float64, g.npts)
    PNL::Array{Complex{T},1} = zeros(Complex{Float64}, Int(g.npts/2+1))
    dEt::Array{T,1} = zeros(Float64, g.npts)
    Ew_0::Array{Complex{T},1} = zeros(Complex{Float64}, Int(g.npts/2+1))
    Ew2::Array{Complex{T},1} = zeros(Complex{Float64}, Int(g.npts/2+1))
    Rt::Array{T,1} = zeros(Float64, g.npts)
    preRt::Array{T,1} = zeros(Float64, g.npts)
    pltype::Symbol
end
