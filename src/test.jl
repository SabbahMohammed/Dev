struct Q
    x::Array{Float64, 1}
end

x = rand(1000);
q = Q(x)

function fx(x::Array{Float64})
    for i in eachindex(x)
        x[i] = exp(-2x[i])
    end
end

function fq(q::Q)
    for i in eachindex(q.x)
        q.x[i] = exp(-2q.x[i])
    end
end

fx1(x) = @time fx(x);
fx1(x);

fq1(q) = @time fq(q);
fq1(q);

