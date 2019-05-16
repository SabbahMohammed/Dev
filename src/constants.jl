#=
constants meaning

τ = the time grid
Et = electric field
Ip = ionization energy for the gas
ns, cn2, wt = constants from adk paper
W = ionization probability
Na = nutral atoms density
Ne = electrons density

p = pressure
Tk = temperture in kelven
c = speed of light
KB = boltezman constant
J = quantum number
h = plancks' constant
e = electron charge
hbar = reduced planck constant
N0 = particles density at atomsperic pressure
m_e = electron mass
ϵ0 = vaccume premitivity
nJ = number of quantum states to be used in Raman QSM
pre = Raman normalization factor obtaied from getpre() function
Bk =The rotational constants
τrot_k = The relaxation factors
=#


const e = 1.60217662e-19
const hbar = 1.0545718e-34
const N0 =  2.5000013629442974e25 # number of atoms in air
const m_e = 9.10938356e-31
const ϵ0 = 8.854e-12
const Tk = 293.0
const c = 3e8
const KB = 1.38064852e-23
const h = 6.62607004e-34

const p_0 = 1 #standered pressure
const Tk_0 = 272.15 #standered temperture
const u = 2.40482555769577 #first order bessel solution http://wwwal.kuicr.kyoto-u.ac.jp/www/accelerator/a4/besselroot.htmlx
const eh = 4.35974434e-18      # Hartree energy
const a_b = 0.52917721092e-10  # Bohr radius
const a_E = 5.14220652e11      # atomic unit of electric field
const a_t = 2.418884326502e-17 # atomic unit of time
