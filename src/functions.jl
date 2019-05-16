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

function get_env(ef)
    ew = fft(ef)
    ew[Int(size(ew)[1]/2+1):end] .= 0.0
    ep = 2.0.*ifft(ew)
    ep
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

function smoothstep(left, right, loss, x)
    x = clamp((x-left)/(right-left), 0, 1)
    return x*x*(3-2*x)*loss
end

function flipsmoothstep(left, right, loss,x)
    if x<left
        return 1.0
    elseif x>right
        return loss
    else
        x = clamp1(1-(x-left)/(right-left), 0, 1.0, loss)
        y = x*x*(3.0-2.0*x)
        return (y*(1-loss)+loss)
    end
end

function clamp1(x, lo, hi, loss)
    if x<lo
        return x = 1.0
    elseif x>hi
        return x = loss
    end
    return x
end


function ml_0!(X, Y, Z, A, z, q)
  k = size(X)[1]
  # @inbounds @simd for j in 1:900
    @inbounds @simd for i in eachindex(X)
      # Ew0[i] = Ew0[i]+1
      X[i] = Y[i]*A[i]*exp(Z[i]*z)*q
    # end
  end
  nothing
end


function ml!(X, Y, z)
  # @inbounds @simd for j in 1:900
    @inbounds @simd for i in eachindex(X)
      # x = X[i]
      X[i] = X[i]*exp(-Y[i]*z)
    # end
  end
  nothing
end

function ml_1!(X, Y, z)
  # @inbounds @simd for j in 1:900
    @inbounds @simd for i in eachindex(X)
      # Ew0[i] = Ew0[i]+1
      X[i] = Y[i]*z
    # end
  end
  nothing
end

function ml_2!(X, Y, Z, z)
  # @inbounds @simd for j in 1:900
    @inbounds @simd for i in eachindex(X)
      # Ew0[i] = Ew0[i]+1
      X[i] = Y[i]*Z[i]*z
    # end
  end
  nothing
end




function ift!(Et, ift, convw)
  # @inbounds @simd for j in 1:900
    @inbounds @simd for i in eachindex(Et)
      # Ew0[i] = Ew0[i]+1
      Et[i] = ift[i]*convw
    # end
  end
  nothing
end

δ(k::Integer,j::Integer) = k == j ? 1 : 0


function cumtrapz(x::AbstractVector, y::AbstractVector, retarr::AbstractVector)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    for i in 2 : length(y)
        retarr[i] = retarr[i-1] + (x[i] - x[i-1]) * (y[i] + y[i-1])
    end
    for i in eachindex(retarr)
        retarr[i] = retarr[i]*1//2
    end
end

function trapz(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y) "x and y vectors must be of the same length!"
    retval = 0.0
    for i in 1 : length(y)-1
        retval += (x[i+1] - x[i]) * (y[i] + y[i+1])
    end
    return 1//2 * retval
end
