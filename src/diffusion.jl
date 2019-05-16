using Plots


function diffusion(t)
    x = range(-0.12,0.12,length = 2^14)
    Ayz = pi*(13e-6)^2
    M = 1e23
    D = 2.14e-6
    return M/(Ayz*sqrt.(4pi*D.*t)).*exp.(-x.^2.0./(4*D.*t))
end

plot(diffusion(50))
#
# function test(t)
#     x = range(-0.12,0.12,length = 2^14)
#     dx = x[2] - x[1]
#     k = 2pi.*fftshift(fftfreq(length(x), 1/dx))
#     Ayz = pi*(13e-6)^2
#     M = 1e23
#     t0 = 0.1
#     D = 2.14e-6
#     C = 0.04
#     u0 = diffsion0(x,t0, M, D, Ayz, 0)
#     uw0 = fftshift(fft(u0))
#     uw = @.  uw0*exp(-k^2*D*t)
#     ux = real.(ifft(ifftshift(uw)))
#     ux
#     prob = ODEProblem(func, uw0, ts, (D, Progress(floor(Int, Float64(0.24*1000)),2)))
#     sol = solve(prob, Vern7(), reltol=1e-4, abstol=1e-10)
#     res = zeros(Complex{Float64}, size(sol))
#     for i in eachindex(sol.t)
#         res[:,i] = ifft(sol[:,i])
#     end
#     res
# end
