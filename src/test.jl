# loss = zeros(size(g.W))
# if fiber.loss == true
#     for i in eachindex(gases2)
#         if gases2[i].type == :O3
#             path = joinpath(pwd(), "ozone")
#             data = readdlm("$path\\full_absorption_data.txt")
#             data = data[:,1]
#             freq = readdlm("$path\\freqdata_ozone.txt")
#             freq = freq[:,1]
            
#             λ = 2pi*3e8./freq*1e9
#             λ_bound = 2.0 .< λ .< 400.0 
#             data = data[λ_bound]
#             freq = freq[λ_bound]
#             data *= gases2[i].p*gases2[i].pp*N0
#             spl = Spline1D(freq, data)
#             loss .= spl.(g.W)
#         end
#     end
# end

# λ = 2pi*3e8./g.W*1e9
# λ_bound = 2.0 .< λ .< 400.0 


# open(dir*"\\loss.txt", "w") do io
#     writedlm(io,  loss)
# end
