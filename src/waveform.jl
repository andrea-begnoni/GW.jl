module waveform
### This is the first module in the pipeline (after the catalog generation), it computes the waveform

# The waveforms presented here are adapted and modified from LALSimulation and GWFAST (https://github.com/CosmoStatGW/gwfast)


import ..UtilsAndConstants as uc    

### Import Julia packages relevant for this module
using DelimitedFiles
using Interpolations
using ForwardDiff
using Roots

export TaylorF2, PhenomD, PhenomD_NRTidal, PhenomHM, PhenomNSBH, PhenomXAS, Model

export Ampl, Phi, Pol, _npar, _event_type, _available_waveforms, _fcut, _finalspin, _radiatednrg, _tau_star, hphc

# Define an abstract type for the models
abstract type Model end

##############################################################################
#   Define general methods a waveform should have
##############################################################################

# Define a function to give an error message if the model is not implemented
function Ampl(
    model::Model,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc
)
    # Implementation specific to each model
    error("Ampl not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function Phi(
    model::Model,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = 0.2,
    fInsJoin_PHI = 0.018,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc
)
    # Implementation specific to each model
    error("Phi not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function Pol(
        model::Model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL
)
    # Implementation specific to each model
    error("Phi not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function _npar(model::Model)
    # Implementation specific to each model
    error("_npar not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define a function to give an error message if the model is not implemented
function _event_type(model::Model) 
   # Implementation specific to each model
   error("_event_type not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

##############################################################################
#   available waveforms
##############################################################################

# Define concrete types for each model
struct PhenomD <: Model end

struct PhenomHM <: Model end

struct TaylorF2 <: Model end

struct PhenomD_NRTidal <: Model end

struct PhenomNSBH <: Model end

struct PhenomXAS <: Model end

function _available_waveforms()
    return ["TaylorF2", "PhenomD", "PhenomHM", "PhenomD_NRTidal", "PhenomNSBH", "PhenomXAS"]
end


#@doc "Function to check the available waveforms and return the corresponding model."
"""
Function to check the available waveforms and return the corresponding model.
"""
function _available_waveforms(waveform::String)

    if waveform == "TaylorF2"
        return TaylorF2()
    elseif waveform == "PhenomD"
        return PhenomD()
    elseif waveform == "PhenomHM"
        return PhenomHM()
    elseif waveform == "PhenomD_NRTidal"
        return PhenomD_NRTidal()
    elseif waveform == "PhenomNSBH"
        return PhenomNSBH()
    elseif waveform == "PhenomXAS"
        return PhenomXAS()
    else
        error("Waveform not available. Choose between: TaylorF2, PhenomD, PhenomHM, PhenomD_NRTidal, PhenomNSBH, PhenomXAS")
    end
end

##############################################################################
#   STRUCTURE USED IN THE MODULE
##############################################################################

struct TF2EccCoeffsStruct
    zero::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    twoV::Union{Float64,ForwardDiff.Dual}
    twoV0::Union{Float64,ForwardDiff.Dual}
    threeV::Union{Float64,ForwardDiff.Dual}
    threeV0::Union{Float64,ForwardDiff.Dual}
    fourV4::Union{Float64,ForwardDiff.Dual}
    fourV2V02::Union{Float64,ForwardDiff.Dual}
    fourV04::Union{Float64,ForwardDiff.Dual}
    fiveV5::Union{Float64,ForwardDiff.Dual}
    fiveV3V02::Union{Float64,ForwardDiff.Dual}
    fiveV2V03::Union{Float64,ForwardDiff.Dual}
    fiveV05::Union{Float64,ForwardDiff.Dual}
    sixV6::Union{Float64,ForwardDiff.Dual}
    sixV4V02::Union{Float64,ForwardDiff.Dual}
    sixV3V03::Union{Float64,ForwardDiff.Dual}
    sixV2V04::Union{Float64,ForwardDiff.Dual}
    sixV06::Union{Float64,ForwardDiff.Dual}
end

struct TF2coeffsStructure
    zero::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
    three::Union{Float64,ForwardDiff.Dual}
    four::Union{Float64,ForwardDiff.Dual}
    five::Union{Float64,ForwardDiff.Dual}
    five_log::Union{Float64,ForwardDiff.Dual}
    six::Union{Float64,ForwardDiff.Dual}
    six_log::Union{Float64,ForwardDiff.Dual}
    seven::Union{Float64,ForwardDiff.Dual}
end

struct PhiInspcoeffsStructure
    initial_phasing::Union{Float64,ForwardDiff.Dual}
    two_thirds::Union{Float64,ForwardDiff.Dual}
    third::Union{Float64,ForwardDiff.Dual}
    third_log::Union{Float64,ForwardDiff.Dual}
    log::Union{Float64,ForwardDiff.Dual}
    min_third::Union{Float64,ForwardDiff.Dual}
    min_two_thirds::Union{Float64,ForwardDiff.Dual}
    min_one::Union{Float64,ForwardDiff.Dual}
    min_four_thirds::Union{Float64,ForwardDiff.Dual}
    min_five_thirds::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    four_thirds::Union{Float64,ForwardDiff.Dual}
    five_thirds::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
end


struct AcoeffsStructure
    two_thirds::Union{Float64,ForwardDiff.Dual}
    one::Union{Float64,ForwardDiff.Dual}
    four_thirds::Union{Float64,ForwardDiff.Dual}
    five_thirds::Union{Float64,ForwardDiff.Dual}
    two::Union{Float64,ForwardDiff.Dual}
    seven_thirds::Union{Float64,ForwardDiff.Dual}
    eight_thirds::Union{Float64,ForwardDiff.Dual}
    three::Union{Float64,ForwardDiff.Dual}
end

struct PhislmpStructure
    two_one::Union{Float64, Vector{Float64}}
    two_two::Union{Float64, Vector{Float64}}
    three_two::Union{Float64, Vector{Float64}}
    three_three::Union{Float64, Vector{Float64}}
    four_three::Union{Float64, Vector{Float64}}
    four_four::Union{Float64, Vector{Float64}}
end

struct AmplslmpStructure
    two_one::Union{Float64, Vector{Float64}}
    two_two::Union{Float64, Vector{Float64}}
    three_two::Union{Float64, Vector{Float64}}
    three_three::Union{Float64, Vector{Float64}}
    four_three::Union{Float64, Vector{Float64}}
    four_four::Union{Float64, Vector{Float64}}
end

##############################################################################
#   Inclusion of the various waveforms
##############################################################################

include("waveform_models/taylorf2.jl")
include("waveform_models/phenomd.jl")
include("waveform_models/phenomd_nrtidal.jl")
include("waveform_models/phenomhm.jl")
include("waveform_models/phenomnsbh.jl")
include("waveform_models/phenomxas.jl")

end# of module
