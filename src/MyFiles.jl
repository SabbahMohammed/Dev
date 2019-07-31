module MyFiles

using SpecialFunctions, Dierckx, SymEngine, DelimitedFiles, ProgressMeter, DSP, FFTW, Roots, Parameters, JSON, NumericalIntegration, Dates
import DifferentialEquations: ODEProblem, solve, Tsit5

export  gnlse, Fiber, createTimeGrid, Pulse, gaussianPulse, sechPulse,
        pulseFromFile, FWHM, Pulse, getElectricField, spectrogram, wbeta2,
        lbeta2, lbeta3, sorder, dlength, L_operator, trapz, beta0, beta1,
        beta2, lfiss, χ3, bound, adk, geisslerPlasma, plasmaTerms, gnlsetest,
        wbeta0, wbeta1, uppe, getElectricFieldE, mfft, mifft, lbeta0,
        uppe2, LE_operator, βnl, dis_wave, zero_disp, lbeta2,
        soliton_order, getGamma, Waist, systemInfo,  soliton_period, smoothstep,
        flipsmoothstep, diffsion, test, diffsion0, rhs2, testrhs, ionRateAdk,
        geisslerPlasmaUPPE, cumtrapz, arb_f, chirpSechPulse, System,
        c, ϵ0, rindex, χ30,refractive_index, O3_disp, dispersion_coeffecient,
        gas, Tk_0, grid, Realχ, splbeta, Realχ, testdis, beta

function __init__()
    global x = symbols(:x)
end

include("constants.jl")
include("grid.jl")
include("raman.jl")
include("functions.jl")

include("gas.jl")
include("plasma.jl")
include("pulse.jl")
include("fiber.jl")
include("dispersion.jl")
include("gnlse.jl")
include("system.jl")

include("uppe.jl")
include("test.jl")


end # module
