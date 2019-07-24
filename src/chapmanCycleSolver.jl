using DifferentialEquations, LinearAlgebra
using Plots
gr()

const N0 =  2.5000013629442974e25
const M0 = 1.5264213305629408e16
# number of air molecules inside the fiber is 1.5264213305629408e16
#=
ozone from
http://www.columbia.edu/itc/chemistry/chem-c2407/hw/ozone_kinetics.pdf
equations 5 & 6
=#



function gaussianPulse(;τ= τ, FWHM=FWHM, t0=t0)
    res0 = sqrt.((exp.(-((τ.-t0).^2.0) ./ (FWHM^2.0))))
end

function createTimeGrid(time_window; npts = 2^14)
    dt = time_window/npts
    τ = LinRange(0, time_window, npts)
    τ = collect(τ)
end

t=createTimeGrid(200)
FWHM=10
gauss(τ) = gaussianPulse(τ= τ, FWHM=FWHM, t0=10) .+ gaussianPulse(τ= τ, FWHM=FWHM, t0=60) .+ gaussianPulse(τ= τ, FWHM=FWHM, t0=110) .+
            gaussianPulse(τ= τ, FWHM=FWHM, t0=160)
# plot(t, gauss.(t))



function fun(du,u,p,t)
    # k1,k2,k3,k4 = p
    k1 = 1.03E-14*gauss(t)
    k3 = 8.58E-2*gauss(t)

    # k1 = 1.03E-14
    k2 = 5.60E-34
    # k3 = 8.58E-2
    # k4 = 8.34E-15
    k4=8.34E-19#sice the reactions are slow compared to the time scale of the ultrafast pulses
    # A = N0*0.50*1e-6*1e-3
    M = M0
    # print(t)
    #=
    u[1] = O, u[2] = O3, u[3] = O2, du[4] = dM = 0
    =#
    du[1] = 2*k1*u[3]-k2*u[1]*u[3]*M+k3*u[2]-k4*u[1]*u[2]
    du[2] = k2*u[1]*u[3]*M - k3*u[2]-k4*u[1]*u[2]
    du[3] = -k1*u[3] - k2*u[3]*M*u[1] + k3*u[2] + 2*k4*u[1]*u[2]
end


# k1 = 1.03E-14
# k2 = 5.60E-34
# k3 = 8.58E-2
# k4 = 8.34E-15
# # A = N0*0.50*1e-6*1e-3
# M = N0*0.80*1e-6*1e-3
p = []

tspan = (0.0, 200)
u0 = [1e7,1e7, N0*0.20*1e-6*1e-3]



prob1 = ODEProblem(fun, u0, tspan, p)
sol = solve(prob1, Tsit5(),reltol = 1e-10, abstol=1e-16)


plot(sol,vars=(0,1))
plot!(sol,vars=(0,2))
# plot!(sol,vars=(0,3))
