function shifted_fft(x)
    return fftshift(fft(fftshift(x)))
end

function shifted_ifft(x)
    return fftshift(ifft(fftshift(x)))
end

function trapz(x, y; dim=1)
    perm = [dim:max(ndims(y),dim); 1:dim-1];
    y = permutedims(y, perm);
    if ndims(y) == 1
        m = 1;
    else
        m = size(y,1);
    end

    if m == 1
        M = length(y);
        out = sum(diff(x).*(y[1:M-1] + y[2:M])/2.0);
    else
        out = transpose(diff(x)) * (reshape(y,Val(2))[1:m-1,:] + reshape(y,Val(2))[2:m,:])/2.0;
        siz = size(y); siz = collect(siz); siz[1] = 1;
        out = reshape(out, tuple(siz...));
        out = permutedims(out, invperm(perm));
        ind = findall(collect(size(out)).==1);
        out = dropdims(out,dims=ind[1]);
        if length(out) == 1;
            out = out[1];
        end
    end

    return out
end


function notfinite(x)
    return ~isfinite(x)
end

function f0(q)
    #=
    same as
    W[~sp.isfinite(W)] = 0.0
    example
    WE = f0(W)
    =#
    # for (i,j) in enumerate(q)
    #     if isfinite(j) == false
    #         q[i]::Float64 = 0
    #     end
    # end
    # q
    if ~isfinite(q)
        return 0
    else
        return q
    end
end

function f1(q)
    #=
    same as
    Na[(~sp.isfinite(Na)) | (Na < 0.0)] = 1.0

    example:
    Ne = [f1(i,0) for i in Na]
    =#
    if isfinite(q) == false || q < 0
        q = 1.0
    end
    q
end

function f10(q)
    #=
    same as
    Na[(~sp.isfinite(Na)) | (Na < 0.0)] = 1.0

    example:
    Ne = [f1(i,0) for i in Na]
    =#
    if isfinite(q) == false || q < 0
        q = 0.0
    end
    q
end

function f2(para1, para2)
    #=
    same as RT[τ<0] = 0
    example:

    raman = sdo(τ, 62.5e-15, 77e-15)
    f2(raman, τ)
    plot(τ, raman)

    =#
    for i in 1:length(para1)
        if para2[i] < 0
            para1[i] = 0
        end
    end
    para1
end

function f3(x)
    return 2pi*3e8/x
end

function isPositive(para)
    if para < 0
        return 0
    end
end

function getElectricField(pulse, fiber)
    ω0 = 2pi*c/pulse.λ0
    Et = pulse.data .* exp.(-1im.*ω0 .* pulse.τ)
    Et .= real(Et)
    Et = sqrt(2/(ϵ0*c*effectiveArea(fiber.diameter))).*Et
    return Et
end

function flip(arr)
    arr1 = zeros(Complex{Float64}, (size(arr)[2], size(arr)[1]))
    for i in 1:size(arr)[2]
        arr1[i,:] = arr[:,i]
    end
    arr1
end

function bound(tmin, tmax, x)
    if x>=tmin && x<=tmax return true
    else return false
    end
end

function spectrogram(τ, U; window=1e-15)
    g(t) = exp.(-1.0*((τ.-t)./(2*window)).^2)
    points = floor(Int, (τ[end] - τ[1] )/window)+1
    out = zeros(Complex{Float64}, (length(τ), points*500))
    counter = minimum(τ)
    for i in 1:points*500
        out[:,i] = shifted_ifft((U .* g(counter)))
        counter += window/500
    end
    #out1 = flip(out)
    abs2.(out)
end

function sorder(γ, P0, τ0, β2)
    return sqrt(γ*P0*τ0/abs(β2))
end

function dlength(τ0, β2)
    return τ0^2/abs(β2)
end

function nlength(γ, P0)
    return 1/(γ*P0)
end


function getElectricFieldE(pulse, fiber)
    ω0 = 2pi*c/pulse.λ0
    Et = pulse.data .* cos.(ω0 .* pulse.τ)
    Et = sqrt(2/(ϵ0*c*effectiveArea(fiber.diameter))).*Et
    return Et
end

function mfft(x)
    y = fft(x)
    if isreal(x[1])
        return y[1:Int(length(y)/2 + 1)]
    else
        return y
    end
end

function mifft(x)
    if length(x) % 2 == 0
        x = cat(x,conj(x[end-1:-1:2]),dims=1)
        return real.(ifft(x))
    else
        return ifft(x)
    end
end



#=
function cumtrapz(x, y; dim=1)
    perm = [dim:max(length(size(y)),dim); 1:dim-1];
    y = permutedims(y, perm);
    if ndims(y) == 1
        n = 1;
        m = length(y);
    else
        m, n = size(y);
    end

    if n == 1
        dt = diff(x)/2.0;
        z = [0; cumsum(dt.*(y[1:(m-1)] + y[2:m]))];
    else
        dt = repeat(diff(x)/2.0,1,n);
        z = [zeros(1,n); cumsum(dt.*(y[1:(m-1), :] + y[2:m, :]),dims=1)];
        z = permutedims(z, invperm(perm));
    end

    return z
end

function cumtrapz0(x, y)
    α = cumtrapz(x, y)
    res = vcat([0], α)
    res
end
=#
