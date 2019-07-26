@with_kw struct GasPlasma{T}
    foo::Int = 6
    ns::Float64
    wp::Float64
    cn2::Float64
    npts::Int64
    Pe::Array{T,1} = zeros(Float64, npts)
    Pe0::Array{T,1} = zeros(Float64, npts)
    Na::Array{T,1} = zeros(Float64, npts)
    Ne::Array{T,1} = zeros(Float64, npts)
    WE::Array{T,1} = zeros(Float64, npts)
    wt::Array{T,1} = zeros(Float64, npts)
    Pe_l::Array{T,1} = zeros(Float64, npts)
    Pe_f::Array{T,1} = zeros(Float64, npts)
    temp1::Array{T,1} = zeros(Float64, npts)
    temp2::Array{T,1} = zeros(Float64, npts)
    temp3::Array{T,1} = zeros(Float64, npts)
    temp4::Array{T,1} = zeros(Float64, npts)
end

struct Gas{T, S}
    type::S #gas type
    p::Float64 #total pressure
    Tk::Float64
    pp::Float64 #partial pressure
    Ip::Float64 #ionization energy
    Reχ::Array{T,1} # refractive index^2 -1
    χ3::Float64
    n2::Float64
    gp::GasPlasma
    #raman to be added
end

function gas(type::Symbol, pp::Float64, Tk::Float64, p::Float64, grid::Grid)
    if type == :Ar
        Ip = 15.75962
    elseif type == :Kr
        Ip = 13.99961
    elseif type == :Ne
        Ip = 21.5646
    elseif type == :He
        Ip = 24.58741
    elseif type == :N2
        Ip = 15.581
    elseif type == :Xe
        Ip = 12.1298
    elseif type == :O2
        Ip = 12.0697
    elseif type == :O3
        Ip = 12.43
    end
    Ip = Ip*e

    Reχ = Realχ(type, p, pp, Tk) # not for O2
    χ3 = χ30(type, pp, Tk)*p
    n2 = n_2(χ3)
    gasplasma = gasPlasma(Ip, grid.npts)
    return Gas{Float64, Symbol}(type, p, Tk, pp, Ip, Reχ.(grid.W), χ3, n2, gasplasma)
end
