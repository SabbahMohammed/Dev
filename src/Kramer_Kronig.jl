function kkrebook2(omega, imchi, alpha)
    if size(omega,1)>size(omega,2)
        omega = omega'
    end
    if size(imchi,1)>size(imchi,2)
        imchi=imchi'
    end
    g = size(omega,2)
    rechi = zeros(size(imchi))
    a = zeros(size(imchi))
    b = zeros(size(imchi))
    deltaomega = omega[2]-omega[1]
    j = 1
    beta1 = 0

    for k=2:g
        b[1] = beta1+imchi[k]*omega[k]^(2*alpha+1)/(omega[k]^2-omega[1]^2)
        beta1 = b[1]
    end
    rechi[1] = 2/pi*deltaomega*b[1]*omega[1]^(-2*alpha)

    j=g
    alpha1 = 0
    for k=1:g-1
        a[g] = alpha1+imchi[k]*omega[k]^(2*alpha+1)/(omega[k]^2-omega[g]^2)
        alpha1 = a[g]
    end
    rechi[g] = 2/pi*deltaomega*a[g]*omega[g]^(-2*alpha)

    for j=2:g-1
        alpha1 = 0
        beta1 = 0
        for k=1:j-1
            a[j] = alpha1+imchi[k]*omega[k]^(2*alpha+1)/(omega[k]^2-omega[j]^2)
            alpha1 = a[j]
        end
        for k=j+1:g
            b[j] = beta1+imchi[k]*omega[k]^(2*alpha+1)/(omega[k]^2-omega[j]^2)
            beta1 = b[j]
        end
        rechi[j] = 2/pi*deltaomega*(a[j]+b[j])*omega[j]^(-2*alpha)
        print("\r$j")
    end
    return rechi
end

imchi = readdlm(raw"C:\Users\ms246\Dropbox (Heriot-Watt University Team)\RES_EPS_Lupo\Projects\Mohammed\phd\ozone absorption spectrum\imchi_0.124nm.txt")
omega = readdlm(raw"C:\Users\ms246\Dropbox (Heriot-Watt University Team)\RES_EPS_Lupo\Projects\Mohammed\phd\ozone absorption spectrum\freqdata_0.124nm.txt")
kk = kkrebook2(omega[2:end], imchi, 0)

dir = raw"C:\Users\ms246\Dropbox (Heriot-Watt University Team)\RES_EPS_Lupo\Projects\Mohammed\phd\ozone absorption spectrum\\"

open(dir*"julia_kk.txt", "w") do io
    writedlm(io, kk)
end

open(dir*"freqJulia.txt", "w") do io
    writedlm(io, omega)
end
