
function Φ(z, λ, w0)
    zR = pi*w0^2/λ
    return atan(z/zR)
end

function w(z, λ, w0)
    zR = pi*w0^2/λ
    return w0*sqrt(1+(z/zR)^2)
end

function R(z, λ, w0)
    zR = pi*w0^2/λ
    return z+zR^2/z
end

function q(z, λ, w0)
    zR = pi*w0^2/λ
    return z + zR*1im
end

function ψ(z, p, m, λ, w0)
    zR = pi*w0^2/λ
    N = abs(m) + 2p
    return (N+1)*atan(z/zR)
end

function Aeff(a)
    return 1.5*a^2
end

# gives fiber output beam waist for pcf
function fiber_beam_waist(a)
    return sqrt(1.5/pi)*a
end

function focal(r, ω0, λ0)
    return sp.pi*ω0*0.64*r/λ0
end


#=
these two functions can be used to design lenses(curved mirrors) system
equivelent to one lens (not available)
=#

function d2(f1,f2, d1)
    #=
    d2 is the distance between the first lens and the beam source.
    be careful with the order of the leses(since the equation is originally for focusing not for collimating).
    For example, the source should always be close to f2(i.e source->space->f2->space->f1->collimated beam)
    these two function work the same way for lenses and curved mirrors
    =#
    return f2*(f1-d1)/(f1+f2-d1)
end

function feff(f1,f2,d1)
    #d1 is the distance between the two lenses
    #these two function work the same way for lenses and curved mirrors
    return (f2*f1)/(f1+f2-d1)
end

#=
Example of ABCD matrex usge.
q1 = q(0, 800e-9, 8.64e-5)
abcd = [1 4; 0 1]*[1 0;-2/2 1]*[1 0.625; 0 1]*[1 0;1/0.5 1]*[1 1.5; 0 1]
the propagation is as follow:
- fiber output with waist w0 = 8.64e-5 (q1)
- 1.5m free space
- concave lens(convex mirror) with focal length(raduis) 0.5(1)
- 0.625m free space
- convex lens(concave mirror) with focal length(raduis) -1(-2)
- 4m free space
this system act as a telescope by magnifiying the initial beam to w1
=#
q1 = q(0, 400e-9, 80e-6)
# abcd = [1 15; 0 1]*[1 0; -2/8 1]*[1 4; 0 1]
abcd = [1 0.381; 0 1]*[1 0;-1/0.4 1]*[1 2.233; 0 1]*[1 0;-1/0.75 1]*[1 0.4; 0 1]*[1 0;1/0.5 1]*[1 2.2; 0 1]
q2 = (abcd[1,1]*q1 + abcd[1,2])/(abcd[2,1]*q1+abcd[2,2])
sol = 1/q2
w1 = sqrt(800e-9/(pi*imag(-sol)))
# w0 = w(2, 800e-9, 8.64e-5)

feff(1,-0.5,0.67)
# d2(0.75,-0.5, 0.342)
