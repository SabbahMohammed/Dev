module MyFiles

using SpecialFunctions, Dierckx, SymEngine, DelimitedFiles, DifferentialEquations,
      ProgressMeter, DSP, FFTW, NumericalIntegration, Roots, Parameters, JSON

export  gnlse, Fiber, createTimeGrid, Pulse, gaussianPulse, sechPulse,
        pulseFromFile, FWHM, Pulse, getElectricField, spectrogram, wbeta2,
        lbeta2, lbeta3, sorder, dlength, L_operator, trapz, beta0, beta1,
        beta2, lfiss, χ3, bound, adk, geisslerPlasma, plasmaTerms, gnlsetest,
        wbeta0, wbeta1, uppe, getElectricFieldE, mfft, mifft, lbeta0,
        test, uppe2, LE_operator, βnl, dis_wave, zero_disp, lbeta2,
        soliton_order, getGamma, Waist, systemInfo, c, soliton_period

function __init__()
    global x = symbols(:x)
end

include("constants.jl")
include("raman.jl")
include("functions.jl")
include("plasma.jl")
include("pulse.jl")
include("dispersion.jl")
include("fiber.jl")
include("gnlse.jl")
#include("test.jl")
include("uppe.jl")


end # module
