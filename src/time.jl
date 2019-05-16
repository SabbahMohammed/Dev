x = rand(4)
y = rand(4)
z = rand(4)
a = zeros(4)



function xx!(a,x,y,z)

    @inbounds @simd for i in 1:4
        a[i] = x[i]*y[i]*exp(z[i])
    end
    nothing
end

function qw()
    @time @. a = x*y*exp(z)
end

function qq()
    @time xx!(a,x,y,z)
end



qw()
qq()
