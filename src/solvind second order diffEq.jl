using DifferentialEquations, LinearAlgebra
using Plots
plotly()

#=
m*x`` + b*x` + k*x + a*x^3 = -mg
y = x`
=#

function fun(du,u,p,t)
    m,k,g,a,b = p
    du[1] = u[2]
    du[2] = -1/m*(b*u[2]+ k*u[1]+ a*u[1]^3+ m*g)
end

m, k, g, a = 1220.0, 35600.0, 17.5, 450000.0
b= 1000.0

tspan = (0.0, 10.0)
u0 = [0.0, 5.0]

prob1 = ODEProblem(fun, u0, tspan, (m, k, g, a, b))
sol = solve(prob1, Tsit5(),reltol = 1e-4, abstol=1e-10)


plot(sol)
