module waveform
### This is the first module in the pipeline, it computes the waveform

# every module imports the previous one and the UtilsAndConstants module
# include("utils.jl")
import ..UtilsAndConstants as uc     #Andrea import vs using

### Import Julia packages relevant for this module
using DelimitedFiles
using Interpolations
using ForwardDiff

export TaylorF2, PhenomD, PhenomD_NRTidal, PhenomHM, PhenomNSBH, Model

export Ampl, Phi, _available_waveforms, _fcut, _finalspin, _radiatednrg, _tau_star, hphc
# Define an abstract type for the models
abstract type Model end


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
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
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
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)
    # Implementation specific to each model
    error("Phi not implemented for model: $(typeof(model)), or there is an error with the number of input parameters")
end

# Define concrete types for each model
struct PhenomD <: Model end

struct PhenomHM <: Model end

struct TaylorF2 <: Model end

struct PhenomD_NRTidal <: Model end

struct PhenomNSBH <: Model end

function _available_waveforms()
    return ["TaylorF2", "PhenomD", "PhenomHM", "PhenomD_NRTidal", "PhenomNSBH"]
end


@doc "Function to check the available waveforms and return the corresponding model."
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
    else
        error("Waveform not available. Choose between: TaylorF2, PhenomD, PhenomHM, PhenomD_NRTidal, PhenomNSBH")
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
# TAYLORF2 3.5 RESTRICTED PN WAVEFORM
##############################################################################


"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(TaylorF2(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The phase is given in radians.

#### Example:
```julia
    f = 10 .^range(1, stop=3, length=100)
    Phi(TaylorF2(), f, 30., 0.25, 0.5, 0.5)
```
"""
function Phi(model::TaylorF2,
    f,
    mc,
    eta,
    chi1,
    chi2,
    Lambda1=0.,
    Lambda2=0.;
    GMsun_over_c3 = uc.GMsun_over_c3,
    is_tidal = false,
    use_QuadMonTid = false,
    is_eccentric = false,
    phiref_vlso = false,
    use_3p5PN_SpinHO = false,
    fRef_ecc = nothing,
)

    # From A. Buonanno, B. Iyer, E. Ochsner, Y. Pan, B.S. Sathyaprakash - arXiv:0907.0700 - eq. (3.18) plus spins as in arXiv:1107.1267 eq. (5.3) up to 2.5PN and PhysRevD.93.084054 eq. (6) for 3PN and 3.5PN
    Mtot_sec = mc * GMsun_over_c3/eta^(0.6)
    v = @. (pi*Mtot_sec*f)^(1. /3.)
    eta2 = eta*eta

    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)        
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    
    chi12, chi22 = chi1*chi1, chi2*chi2
    chi1dotchi2  = chi1*chi2
    chi_s, chi_a   = 0.5*(chi1 + chi2), 0.5*(chi1 - chi2)
    chi_s2, chi_a2 = chi_s*chi_s, chi_a*chi_a
    chi_sdotchi_a  = chi_s*chi_a

    vlso = 1. /sqrt(6.)
    
    if is_tidal && use_QuadMonTid
        # A non-zero tidal deformability induces a quadrupole moment (for BBH it is 1).
        # The relation between the two is given in arxiv:1608.02582 eq. (15) with coefficients from third row of Table I
        # We also extend the range to 0 <= Lam < 1, as done in LALSimulation in LALSimUniversalRelations.c line 123
        QuadMon1 = ifelse(Lambda1 < 1., 1. + Lambda1*(0.427688866723244 + Lambda1*(-0.324336526985068 + Lambda1*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda1) + 0.04812 * log(Lambda1)^2 -4.283e-3 * log(Lambda1)^3 + 1.245e-4 * log(Lambda1)^4))
        QuadMon2 = ifelse(Lambda2 < 1., 1. + Lambda2*(0.427688866723244 + Lambda2*(-0.324336526985068 + Lambda2*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda2) + 0.04812 * log(Lambda2)^2 -4.283e-3 * log(Lambda2)^3 + 1.245e-4 * log(Lambda2)^4))
    else
        QuadMon1, QuadMon2 = 1.,1.
    end
    TF2OverallAmpl = 3. /(128. * eta)
    TF2_5coeff_tmp = (732985. /2268. - 24260. *eta/81. - 340. *eta2/9.)*chi_s + (732985. /2268. + 140. *eta/9.)*Seta*chi_a

    if phiref_vlso
        TF2coeffs_five = (38645. *pi/756. - 65. *pi*eta/9. - TF2_5coeff_tmp)*(1. -3. *log(vlso))
        phiR = 0.
    else
        TF2coeffs_five = (38645. *pi/756. - 65. *pi*eta/9. - TF2_5coeff_tmp)
        # This pi factor is needed to include LAL fRef rescaling, so to end up with the exact same waveform
        phiR = pi
    end

    if use_3p5PN_SpinHO
        # This part includes SS and SSS contributions at 3.5PN, which are not included in LAL
            TF2coeffs_seven = 77096675. *pi/254016. + 378515. *pi*eta/1512. - 74045. *pi*eta2/756. + (-25150083775. /3048192. + 10566655595. *eta/762048. - 1042165. *eta2/3024. + 5345. *eta^3/36. + (14585. /8. - 7270. *eta + 80. *eta2)*chi_a2)*chi_s + (14585. /24. - 475. *eta/6. + 100. *eta2/3.)*chi_s2*chi_s + Seta*((-25150083775. /3048192. + 26804935. *eta/6048. - 1985. *eta2/48.)*chi_a + (14585. /24. - 2380. *eta)*chi_a2*chi_a + (14585. /8. - 215. *eta/2.)*chi_a*chi_s2)
        else
            TF2coeffs_seven = 77096675. *pi/254016. + 378515. *pi*eta/1512. - 74045. *pi*eta2/756. + (-25150083775. /3048192. + 10566655595. *eta/762048. - 1042165. *eta2/3024. + 5345. *eta^3/36.)*chi_s + Seta*((-25150083775. /3048192. + 26804935. *eta/6048. - 1985. *eta2/48.)*chi_a)
        end

    TF2coeffs = TF2coeffsStructure(
        1.,
        0.,
        3715. /756. + (55. *eta)/9.,
        -16. *pi + (113. *Seta*chi_a)/3. + (113. /3. - (76. *eta)/3.)*chi_s,
        5. *(3058.673/7.056 + 5429. /7. *eta+617. *eta2)/72. + 247. /4.8*eta*chi1dotchi2 -721. /4.8*eta*chi1dotchi2 + (-720. /9.6*QuadMon1 + 1. /9.6)*m1ByM^2*chi12 + (-720. /9.6*QuadMon2 + 1. /9.6)*m2ByM^2*chi22 + (240. /9.6*QuadMon1 - 7. /9.6)*m1ByM^2*chi12 + (240. /9.6*QuadMon2 - 7. /9.6)*m2ByM^2*chi22,
        TF2coeffs_five,
        (38645. *pi/756. - 65. *pi*eta/9. - TF2_5coeff_tmp)*3,
        11583.231236531/4.694215680 - 640. /3. *pi^2 - 684.8/2.1* MathConstants.eulergamma + eta*(-15737.765635/3.048192 + 225.5/1.2*pi^2) + eta2*76.055/1.728 - eta^3*127.825/1.296 - log(4.)*684.8/2.1 + pi*chi1*m1ByM*(1490. /3. + m1ByM*260.) + pi*chi2*m2ByM*(1490. /3. + m2ByM*260.) + (326.75/1.12 + 557.5/1.8*eta)*eta*chi1dotchi2 + (4703.5/8.4+2935. /6. *m1ByM-120. *m1ByM^2)*m1ByM^2*QuadMon1*chi12 + (-4108.25/6.72-108.5/1.2*m1ByM+125.5/3.6*m1ByM^2)*m1ByM^2*chi12 + (4703.5/8.4+2935. /6. *m2ByM-120. *m2ByM^2)*m2ByM^2*QuadMon2*chi22 + (-4108.25/6.72-108.5/1.2*m2ByM+125.5/3.6*m2ByM^2)*m2ByM^2*chi22,
        -(6848. /21.),
        TF2coeffs_seven,
    )

    if is_eccentric
        # These are the eccentricity dependent coefficients up to 3 PN order, in the low-eccentricity limit, from arXiv:1605.00304
        
        if isnothing(fRef_ecc)
            v0ecc = amin(v, axis=0)
        else
            v0ecc = (pi*Mtot_sec*fRef_ecc)^(1. /3.)
        end  
        TF2EccCoeffs = TF2EccCoeffsStruct(
            1.,
            0.,
            29.9076223/8.1976608 + 18.766963/2.927736*eta,
            2.833/1.008 - 19.7/3.6*eta,
            -28.19123/2.82600*pi,
            37.7/7.2*pi,
            16.237683263/3.330429696 + 241.33060753/9.71375328*eta+156.2608261/6.9383952*eta2,
            84.7282939759/8.2632420864-7.18901219/3.68894736*eta-36.97091711/1.05398496*eta2,
            -1.193251/3.048192 - 66.317/9.072*eta +18.155/1.296*eta2,
            -28.31492681/1.18395270*pi - 115.52066831/2.70617760*pi*eta,
            -79.86575459/2.84860800*pi + 55.5367231/1.0173600*pi*eta,
            112.751736071/5.902315776*pi + 70.75145051/2.10796992*pi*eta,
            76.4881/9.0720*pi - 94.9457/2.2680*pi*eta,
            -436.03153867072577087/1.32658535116800000 + 53.6803271/1.9782000*MathConstants.eulergamma + 157.22503703/3.25555200*pi^2 +(2991.72861614477/6.89135247360 - 15.075413/1.446912*pi^2)*eta +345.5209264991/4.1019955200*eta2 + 506.12671711/8.78999040*eta^3 + 384.3505163/5.9346000*log(2.) - 112.1397129/1.7584000*log(3.),
            46.001356684079/3.357073133568 + 253.471410141755/5.874877983744*eta - 169.3852244423/2.3313007872*eta2 - 307.833827417/2.497822272*eta^3,
            -106.2809371/2.0347200*pi^2,
            -3.56873002170973/2.49880440692736 - 260.399751935005/8.924301453312*eta + 15.0484695827/3.5413894656*eta2 + 340.714213265/3.794345856*eta^3,
            265.31900578691/1.68991764480 - 33.17/1.26*MathConstants.eulergamma + 12.2833/1.0368*pi^2 + (91.55185261/5.48674560 - 3.977/1.152*pi^2)*eta - 5.732473/1.306368*eta2 - 30.90307/1.39968*eta^3 + 87.419/1.890*log(2.) - 260.01/5.60*log(3.),
        )
        
        TF2EccOverallAmpl = -2.355/1.462*ecc^2*((v0ecc/v)^(19. /3.))
        phi_Ecc = TF2EccOverallAmpl * (TF2EccCoeffs.zero + TF2EccCoeffs.one*v + (TF2EccCoeffs.twoV*v^2 + TF2EccCoeffs.twoV0*v0ecc^2) + (TF2EccCoeffs.threeV*v^3 + TF2EccCoeffs.threeV0*v0ecc^3) + (TF2EccCoeffs.fourV4*v^4 + TF2EccCoeffs.fourV2V02*v^2*v0ecc^2 + TF2EccCoeffs.fourV04*v0ecc^4) + (TF2EccCoeffs.fiveV5*v^5 + TF2EccCoeffs.fiveV3V02*v^3*v0ecc^2 + TF2EccCoeffs.fiveV2V03*v^2*v0ecc^3 + TF2EccCoeffs.fiveV05*v0ecc^5) + ((TF2EccCoeffs.sixV6 + 53.6803271/3.9564000*log(16. *v^2))*(v^6)+ TF2EccCoeffs.sixV4V02*v^4*v0ecc^2 + TF2EccCoeffs.sixV3V03*v^3*v0ecc^3 + TF2EccCoeffs.sixV2V04*v^2*v0ecc^4 + (TF2EccCoeffs.sixV06 - 33.17/2.52*log(16. *v0ecc^2))*v0ecc^6))
    else
        phi_Ecc = 0.
    end    
    if is_tidal
        # Add tidal contribution if needed, as in PhysRevD.89.103012
        Lam_t, delLam    = uc.Lamt_delLam_from_Lam12(Lambda1, Lambda2, eta)
        
        phi_Tidal = (-0.5*39. *Lam_t)*(v^10.) + (-3115. /64. *Lam_t + 6595. /364. *Seta*delLam)*(v^12.)
        
    else
        phi_Tidal = 0.
    end
    phase = @. TF2OverallAmpl*(TF2coeffs.zero + TF2coeffs.one*v + TF2coeffs.two*v^2 + TF2coeffs.three*v^3 + TF2coeffs.four*v^4 + (TF2coeffs.five + TF2coeffs.five_log*log(v))*v^5 + (TF2coeffs.six + TF2coeffs.six_log*log(v))*v^6 + TF2coeffs.seven*v^7 + phi_Tidal + phi_Ecc)/(v^5.)
    
    return phase .+ phiR .- pi*0.25
end

function Ampl(model::TaylorF2, f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)
    return Ampl(model, f, mc, dL, clightGpc = clightGpc, GMsun_over_c3 = GMsun_over_c3)
end
"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(TaylorF2(), f, mc, dL)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `dL`: Luminosity distance to the source, in Gpc.

#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The amplitude is given as a strain and is dimensionless.

#### Example:
```julia
    f = 10 .^range(1, stop=3, length=100)
    Ampl(TaylorF2(), f, 30., 8.)
```
"""
function Ampl(model::TaylorF2, f, mc, dL; clightGpc = uc.clightGpc, GMsun_over_c3 = uc.GMsun_over_c3)

    # In the restricted PN approach the amplitude is the same as for the Newtonian approximation, so this term is equivalent
    amplitude = @. sqrt(5. /24.) * (pi^(-2. /3.)) * clightGpc/dL * (GMsun_over_c3*mc)^(5. /6.) * (f^(-7. /6.))
    return amplitude
end


"""
Compute the cut frequency of the waveform as a function of the events parameters, in Hz.

This can be approximated as 2 f_ISCO for inspiral only waveforms:

    - if Schwarzschild ISCO for a non-rotating final BH is used (default), use _fcut(model::TaylorF2, mc, eta) to compute the cut frequency of the waveform for the chosen events, in Hz.
    - if Kerr ISCO for a rotating final BH is computed, use _fcut(model::TaylorF2, mc, eta, chi1, chi2) to compute the cut frequency of the waveform for the chosen events, in Hz. As in `arXiv:2108.05861 <https://arxiv.org/abs/2108.05861>`_ (see in particular App. C). NOTE: this is pushing the validity of the model to the limit, and is not the default option.
#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `model`: Model of the waveform.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Return:
-  Cut frequency (float) of the waveform for the chosen event, in Hz.

"""
function _fcut(model::TaylorF2, mc, eta, chi1, chi2; GMsun_over_c3 = uc.GMsun_over_c3)

    println("Using fcut of Kerr orbit since the two angular momenta 'chi1' and 'chi2' were given")
    # println("chi1: ", chi1)
    # println("chi2: ", chi2)
    eta2 = eta*eta
    Mtot = mc/(eta^(3. /5.))
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    s = (m1^2 * chi1 + m2^2 * chi2) / (m1^2 + m2^2)
    atot = (chi1 + chi2*(m2/m1)*(m2/m1))/((1. +m2/m1)*(1. +m2/m1))
    aeff = atot + 0.41616*eta*(chi1 + chi2)

    r_ISCO = _r_ISCO_of_chi(aeff)
    
    EradNS = eta * (0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 + 3.352411249771192 * eta^3)
    EradTot = (EradNS * (1. + (-0.0030302335878845507 - 2.0066110851351073 * eta + 7.7050567802399215 * eta2) * s)) / (1. + (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s)
    
    Mfin = Mtot*(1. -EradTot)
    L_ISCO = 2. /(3. *sqrt(3.))*(1. + 2. *sqrt(3. *r_ISCO - 2.))
    E_ISCO = sqrt(1. - 2. /(3. *r_ISCO))
    
    chif = atot + eta*(L_ISCO - 2. *atot*(E_ISCO - 1.)) + (-3.821158961 - 1.2019*aeff - 1.20764*aeff^2)*eta2 + (3.79245 + 1.18385*aeff + 4.90494*aeff^3)*eta^3
    
    Om_ISCO = 1. /(((_r_ISCO_of_chi(chif))^(3. /2.))+chif)

    return Om_ISCO/(pi*Mfin*GMsun_over_c3)
end

function _r_ISCO_of_chi(chi)
    Z1_ISCO = 1.0 + ((1.0 - chi^2)^(1. /3.))*((1.0+chi)^(1. /3.) + (1.0-chi)^(1. /3.))
    Z2_ISCO = sqrt(3.0*chi*chi + Z1_ISCO^2)
    return ifelse(chi>0., 3.0 + Z2_ISCO - sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)), 3.0 + Z2_ISCO + sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)))
end

##############################################################################
# IMRPhenomD WAVEFORM
##############################################################################

@doc raw"""
IMRPhenomD waveform model.

Relevant references:
    [1] `arXiv:1508.07250 <https://arxiv.org/abs/1508.07250>`_
    
    [2] `arXiv:1508.07253 <https://arxiv.org/abs/1508.07253>`_
    All is taken from LALSimulation and arXiv:1508.07250, arXiv:1508.07253

"""



function _readQNMgrid_a(pathGWFASTWF::String)
    return readdlm(pathGWFASTWF * "QNMData_a.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
end

function _readQNMgrid_fring(pathGWFASTWF::String)
    return readdlm(pathGWFASTWF * "QNMData_fring.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
end

function _readQNMgrid_fdamp(pathGWFASTWF::String)
    return readdlm(pathGWFASTWF * "QNMData_fdamp.txt")[:, 1]   # [:,1] is to make it a 1D array instead of a 2D array
end

""" helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Phi(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2,
    Lambda1,
    Lambda2,
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    container = nothing,
)
    # raise error if Lambda1 or Lambda2 are different from zero
    # if Lambda1 != 0.0 || Lambda2 != 0.0
    #     println("You requested Lambda1 and/or Lambda2 different from zero")
    #     println("Tidal parameters Lambda1 and Lambda2 are not supported in PhenomD")
    #     println("The code put them to zero")
    # end
    return Phi(model, f, mc, eta, chi1, chi2, fInsJoin_PHI=fInsJoin_PHI, fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3, container=container)
end

## # Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude
# fInsJoin_Ampl = 0.014
# # Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase
# fInsJoin_PHI = 0.018
# # Dimensionless frequency (Mf) at which we define the end of the waveform
# fcutPar = 0.2
"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomD(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Optional arguments:
-  `fInsJoin_PHI`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The phase is given in radians.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomD(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Phi(PhenomD(), f, mc, eta, chi1, chi2)
```
"""
function Phi(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    container = nothing,
)

    path=pwd()*"/useful_files/WFfiles/"
    QNMgrid_a = _readQNMgrid_a(path)
    QNMgrid_fring = _readQNMgrid_fring(path)
    QNMgrid_fdamp = _readQNMgrid_fdamp(path)

    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta

    pi2 = pi * pi

    fInsJoin = fInsJoin_PHI
    QuadMon1, QuadMon2 = 1.0, 1.0

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)  
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)

    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    #fgrid = M*GMsun_over_c3*f
    fgrid = M * GMsun_over_c3 .* f # ANDREA: maybe it would be more useful if this function works with arrays of sources of size n_processes
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    # Compute ringdown and damping frequencies from interpolators
    fring = LinearInterpolation(QNMgrid_a, QNMgrid_fring)(aeff) / (1.0 - Erad)
    fdamp = LinearInterpolation(QNMgrid_a, QNMgrid_fdamp)(aeff) / (1.0 - Erad)

    # Compute sigma coefficients appearing in arXiv:1508.07253 eq. (28)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    sigma1 =
        2096.551999295543 +
        1463.7493168261553 * eta +
        (
            1312.5493286098522 + 18307.330017082117 * eta - 43534.1440746107 * eta2 +
            (-833.2889543511114 + 32047.31997183187 * eta - 108609.45037520859 * eta2) *
            xi +
            (452.25136398112204 + 8353.439546391714 * eta - 44531.3250037322 * eta2) * xi2
        ) * xi
    sigma2 =
        -10114.056472621156 - 44631.01109458185 * eta +
        (
            -6541.308761668722 - 266959.23419307504 * eta +
            686328.3229317984 * eta2 +
            (3405.6372187679685 - 437507.7208209015 * eta + 1.6318171307344697e6 * eta2) *
            xi +
            (-7462.648563007646 - 114585.25177153319 * eta + 674402.4689098676 * eta2) * xi2
        ) * xi
    sigma3 =
        22933.658273436497 +
        230960.00814979506 * eta +
        (
            14961.083974183695 + 1.1940181342318142e6 * eta - 3.1042239693052764e6 * eta2 +
            (-3038.166617199259 + 1.8720322849093592e6 * eta - 7.309145012085539e6 * eta2) *
            xi +
            (42738.22871475411 + 467502.018616601 * eta - 3.064853498512499e6 * eta2) * xi2
        ) * xi
    sigma4 =
        -14621.71522218357 - 377812.8579387104 * eta +
        (
            -9608.682631509726 - 1.7108925257214056e6 * eta +
            4.332924601416521e6 * eta2 +
            (
                -22366.683262266528 - 2.5019716386377467e6 * eta +
                1.0274495902259542e7 * eta2
            ) * xi +
            (-85360.30079034246 - 570025.3441737515 * eta + 4.396844346849777e6 * eta2) *
            xi2
        ) * xi

    # Compute beta coefficients appearing in arXiv:1508.07253 eq. (16)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    beta1 =
        97.89747327985583 - 42.659730877489224 * eta +
        (
            153.48421037904913 - 1417.0620760768954 * eta +
            2752.8614143665027 * eta2 +
            (138.7406469558649 - 1433.6585075135881 * eta + 2857.7418952430758 * eta2) *
            xi +
            (41.025109467376126 - 423.680737974639 * eta + 850.3594335657173 * eta2) * xi2
        ) * xi
    beta2 =
        -3.282701958759534 - 9.051384468245866 * eta +
        (
            -12.415449742258042 + 55.4716447709787 * eta - 106.05109938966335 * eta2 +
            (-11.953044553690658 + 76.80704618365418 * eta - 155.33172948098394 * eta2) *
            xi +
            (-3.4129261592393263 + 25.572377569952536 * eta - 54.408036707740465 * eta2) *
            xi2
        ) * xi
    beta3 =
        -0.000025156429818799565 +
        0.000019750256942201327 * eta +
        (
            -0.000018370671469295915 +
            0.000021886317041311973 * eta +
            0.00008250240316860033 * eta2 +
            (
                7.157371250566708e-6 - 0.000055780000112270685 * eta +
                0.00019142082884072178 * eta2
            ) * xi +
            (
                5.447166261464217e-6 - 0.00003220610095021982 * eta +
                0.00007974016714984341 * eta2
            ) * xi2
        ) * xi

    # Compute alpha coefficients appearing in arXiv:1508.07253 eq. (14)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    alpha1 =
        43.31514709695348 +
        638.6332679188081 * eta +
        (
            -32.85768747216059 + 2415.8938269370315 * eta - 5766.875169379177 * eta2 +
            (-61.85459307173841 + 2953.967762459948 * eta - 8986.29057591497 * eta2) * xi +
            (-21.571435779762044 + 981.2158224673428 * eta - 3239.5664895930286 * eta2) *
            xi2
        ) * xi
    alpha2 =
        -0.07020209449091723 - 0.16269798450687084 * eta +
        (
            -0.1872514685185499 + 1.138313650449945 * eta - 2.8334196304430046 * eta2 +
            (-0.17137955686840617 + 1.7197549338119527 * eta - 4.539717148261272 * eta2) *
            xi +
            (-0.049983437357548705 + 0.6062072055948309 * eta - 1.682769616644546 * eta2) *
            xi2
        ) * xi
    alpha3 =
        9.5988072383479 - 397.05438595557433 * eta +
        (
            16.202126189517813 - 1574.8286986717037 * eta +
            3600.3410843831093 * eta2 +
            (27.092429659075467 - 1786.482357315139 * eta + 5152.919378666511 * eta2) * xi +
            (11.175710130033895 - 577.7999423177481 * eta + 1808.730762932043 * eta2) * xi2
        ) * xi
    alpha4 =
        -0.02989487384493607 +
        1.4022106448583738 * eta +
        (
            -0.07356049468633846 +
            0.8337006542278661 * eta +
            0.2240008282397391 * eta2 +
            (-0.055202870001177226 + 0.5667186343606578 * eta + 0.7186931973380503 * eta2) *
            xi +
            (
                -0.015507437354325743 +
                0.15750322779277187 * eta +
                0.21076815715176228 * eta2
            ) * xi2
        ) * xi
    alpha5 =
        0.9974408278363099 - 0.007884449714907203 * eta +
        (
            -0.059046901195591035 + 1.3958712396764088 * eta - 4.516631601676276 * eta2 +
            (-0.05585343136869692 + 1.7516580039343603 * eta - 5.990208965347804 * eta2) *
            xi +
            (-0.017945336522161195 + 0.5965097794825992 * eta - 2.0608879367971804 * eta2) *
            xi2
        ) * xi

    # Compute the TF2 phase coefficients and put them in a dictionary (spin effects are included up to 3.5PN)
    TF2OverallAmpl = 3 / (128.0 * eta)

    # For 3PN coeff we use chi1 and chi2 so to have the quadrupole moment explicitly appearing


    TF2_5coeff_tmp =
        38645.0 * pi / 756.0 - 65.0 * pi * eta / 9.0 - (
            (732985.0 / 2268.0 - 24260.0 * eta / 81.0 - 340.0 * eta2 / 9.0) * chi_s +
            (732985.0 / 2268.0 + 140.0 * eta / 9.0) * Seta * chi_a
        ) #variable to be used later
    TF2_6coeff_tmp =
        11583.231236531 / 4.694215680 - 640.0 / 3.0 * pi2 -
        684.8 / 2.1 * MathConstants.eulergamma +
        eta * (-15737.765635 / 3.048192 + 225.5 / 1.2 * pi2) +
        eta2 * 76.055 / 1.728 - eta2 * eta * 127.825 / 1.296 - log(4.0) * 684.8 / 2.1 +
        pi * chi1 * m1ByM * (1490.0 / 3.0 + m1ByM * 260.0) +
        pi * chi2 * m2ByM * (1490.0 / 3.0 + m2ByM * 260.0) +
        (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) *
        m1ByM^2 *
        QuadMon1 *
        chi12 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 ) *
        m1ByM^2 *
        chi12 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) *
        m2ByM^2 *
        QuadMon2 *
        chi22 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 ) *
        m2ByM^2 *
        chi22

    TF2coeffs = TF2coeffsStructure(
        1.0,
        0.0,
        3715.0 / 756.0 + (55.0 * eta) / 9.0,
        -16.0 * pi +
        (113.0 * Seta * chi_a) / 3.0 +
        (113.0 / 3.0 - (76.0 * eta) / 3.0) * chi_s,
        5.0 * (3058.673 / 7.056 + 5429.0 / 7.0 * eta + 617.0 * eta2) / 72.0 +
        247.0 / 4.8 * eta * chi1dotchi2 - 721.0 / 4.8 * eta * chi1dotchi2 +
        (-720.0 / 9.6 * QuadMon1 + 1.0 / 9.6) * m1ByM^2  * chi12 +
        (-720.0 / 9.6 * QuadMon2 + 1.0 / 9.6) * m2ByM^2  * chi22 +
        (240.0 / 9.6 * QuadMon1 - 7.0 / 9.6) * m1ByM^2  * chi12 +
        (240.0 / 9.6 * QuadMon2 - 7.0 / 9.6) * m2ByM^2  * chi22,
        TF2_5coeff_tmp,
        TF2_5coeff_tmp * 3.0,
        TF2_6coeff_tmp - (
            (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 )
            ) *
            m1ByM *
            m1ByM *
            chi12 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 )
            ) *
            m2ByM *
            m2ByM *
            chi22
        ),
        -6848.0 / 21.0,
        77096675.0 * pi / 254016.0 + 378515.0 * pi * eta / 1512.0 -
        74045.0 * pi * eta2 / 756.0 +
        (
            -25150083775.0 / 3048192.0 + 10566655595.0 * eta / 762048.0 -
            1042165.0 * eta2 / 3024.0 + 5345.0 * eta2 * eta / 36.0
        ) * chi_s +
        Seta * (
            (
                -25150083775.0 / 3048192.0 + 26804935.0 * eta / 6048.0 -
                1985.0 * eta2 / 48.0
            ) * chi_a
        ),
    )


    PhiInspcoeffs = PhiInspcoeffsStructure(
        TF2coeffs.five * TF2OverallAmpl,
        TF2coeffs.seven * TF2OverallAmpl * (pi^(2.0 / 3.0)),
        TF2coeffs.six * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.six_log * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.five_log * TF2OverallAmpl,
        TF2coeffs.four * TF2OverallAmpl * (pi^(-1.0 / 3.0)),
        TF2coeffs.three * TF2OverallAmpl * (pi^(-2.0 / 3.0)),
        TF2coeffs.two * TF2OverallAmpl / pi,
        TF2coeffs.one * TF2OverallAmpl * (pi^(-4.0 / 3.0)),
        TF2coeffs.zero * TF2OverallAmpl * (pi^(-5.0 / 3.0)),
        sigma1,
        sigma2 * 0.75,
        sigma3 * 0.6,
        sigma4 * 0.5,
    )

    #Now compute the coefficients to align the three parts

    fMRDJoin = 0.5 * fring

    # First the Inspiral - Intermediate: we compute C1Int and C2Int coeffs
    # Equations to solve for to get C(1) continuous join
    # PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    # Joining at fInsJoin
    # PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    # PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    # This is the first derivative wrt f of the inspiral phase computed at fInsJoin, first add the PN contribution and then the higher order calibrated terms

    DPhiIns =
        (
            2.0 * TF2coeffs.seven * TF2OverallAmpl * ((pi * fInsJoin)^(7.0 / 3.0)) +
            (
                TF2coeffs.six * TF2OverallAmpl +
                TF2coeffs.six_log * TF2OverallAmpl * (1.0 + log(pi * fInsJoin) / 3.0)
            ) * ((pi * fInsJoin)^(2.0)) +
            TF2coeffs.five_log * TF2OverallAmpl * ((pi * fInsJoin)^(5.0 / 3.0)) -
            TF2coeffs.four * TF2OverallAmpl * ((pi * fInsJoin)^(4.0 / 3.0)) -
            2.0 * TF2coeffs.three * TF2OverallAmpl * (pi * fInsJoin) -
            3.0 * TF2coeffs.two * TF2OverallAmpl * ((pi * fInsJoin)^(2.0 / 3.0)) -
            4.0 * TF2coeffs.one * TF2OverallAmpl * ((pi * fInsJoin)^(1.0 / 3.0)) -
            5.0 * TF2coeffs.zero * TF2OverallAmpl
        ) * pi / (3.0 * ((pi * fInsJoin)^(8.0 / 3.0))) +
        (
            sigma1 +
            sigma2 * (fInsJoin^(1.0 / 3.0)) +
            sigma3 * (fInsJoin^(2.0 / 3.0)) +
            sigma4 * fInsJoin
        ) * etaInv

    # This is the first derivative of the Intermediate phase computed at fInsJoin
    DPhiInt = (beta1 + beta3 / (fInsJoin^4) + beta2 / fInsJoin) * etaInv
    C2Int = DPhiIns - DPhiInt

    # This is the inspiral phase computed at fInsJoin
    # PhiInsJoin = PhiInspcoeffs['initial_phasing'] + PhiInspcoeffs['two_thirds']*(fInsJoin^(2. /3.)) + PhiInspcoeffs['third']*(fInsJoin^(1. /3.)) + PhiInspcoeffs['third_log']*(fInsJoin^(1. /3.))*log(pi*fInsJoin)/3. + PhiInspcoeffs['log']*log(pi*fInsJoin)/3. + PhiInspcoeffs['min_third']*(fInsJoin^(-1. /3.)) + PhiInspcoeffs['min_two_thirds']*(fInsJoin^(-2. /3.)) + PhiInspcoeffs['min_one']/fInsJoin + PhiInspcoeffs['min_four_thirds']*(fInsJoin^(-4. /3.)) + PhiInspcoeffs['min_five_thirds']*(fInsJoin^(-5. /3.)) + (PhiInspcoeffs['one']*fInsJoin + PhiInspcoeffs['four_thirds']*(fInsJoin^(4. /3.)) + PhiInspcoeffs['five_thirds']*(fInsJoin^(5. /3.)) + PhiInspcoeffs['two']*fInsJoin*fInsJoin)/eta
    PhiInsJoin =
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (fInsJoin^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (fInsJoin^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (fInsJoin^(1.0 / 3.0)) * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.log * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.min_third * (fInsJoin^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (fInsJoin^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / fInsJoin +
        PhiInspcoeffs.min_four_thirds * (fInsJoin^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (fInsJoin^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * fInsJoin +
            PhiInspcoeffs.four_thirds * (fInsJoin^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (fInsJoin^(5.0 / 3.0)) +
            PhiInspcoeffs.two * fInsJoin * fInsJoin
        ) * etaInv
    # This is the Intermediate phase computed at fInsJoin
    PhiIntJoin =
        beta1 * fInsJoin - beta3 / (3.0 * fInsJoin * fInsJoin * fInsJoin) +
        beta2 * log(fInsJoin)

    C1Int = PhiInsJoin - PhiIntJoin * etaInv - C2Int * fInsJoin

    # Now the same for Intermediate - Merger-Ringdown: we also need a temporary Intermediate Phase function
    PhiIntTempVal =
        (beta1 * fMRDJoin - beta3 / (3.0 * fMRDJoin^3) + beta2 * log(fMRDJoin)) * etaInv +
        C1Int +
        C2Int * fMRDJoin
    DPhiIntTempVal = C2Int + (beta1 + beta3 / (fMRDJoin^4) + beta2 / fMRDJoin) * etaInv
    DPhiMRDVal =
        (
            alpha1 +
            alpha2 / (fMRDJoin^2) +
            alpha3 / (fMRDJoin^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fMRDJoin - alpha5 * fring) * (fMRDJoin - alpha5 * fring) / (fdamp^2)
                )
            )
        ) * etaInv
    PhiMRJoinTemp =
        -(alpha2 / fMRDJoin) +
        (4.0 / 3.0) * (alpha3 * (fMRDJoin^(0.75))) +
        alpha1 * fMRDJoin +
        alpha4 * atan((fMRDJoin - alpha5 * fring) / fdamp)

    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRD * fMRDJoin
    # Time shift so that peak amplitude is approximately at t=0
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi^2
        ) * xi
    gamma3 =
        1.3081615607036106 - 0.005537729694807678 * eta +
        (
            -0.06782917938621007 - 0.6689834970767117 * eta +
            3.403147966134083 * eta2 +
            (-0.05296577374411866 - 0.9923793203111362 * eta + 4.820681208409587 * eta2) *
            xi +
            (
                -0.006134139870393713 - 0.38429253308696365 * eta +
                1.7561754421985984 * eta2
            ) * xi^2
        ) * xi
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs.(fring - (fdamp * gamma3) / gamma2),
        abs.(fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2^2)) * gamma3) / gamma2),
    )
    t0 =
        (
            alpha1 +
            alpha2 / (fpeak^2) +
            alpha3 / (fpeak^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fpeak - alpha5 * fring) * (fpeak - alpha5 * fring) / (fdamp * fdamp)
                )
            )
        ) * etaInv

    # LAL sets fRef as the minimum frequency, do the same
    fRef = fgrid[1]   #minimum(fgrid)

    phiRef = ifelse(
        fRef < fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fRef^(2.0 / 3.0) +
        PhiInspcoeffs.third * fRef^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fRef^(1.0 / 3.0) * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.log * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.min_third * fRef^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fRef^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fRef +
        PhiInspcoeffs.min_four_thirds * fRef^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fRef^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fRef +
            PhiInspcoeffs.four_thirds * fRef^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fRef^(5.0 / 3.0) +
            PhiInspcoeffs.two * fRef * fRef
        ) * etaInv,
        ifelse(
            fRef < fMRDJoin,
            (beta1 * fRef - beta3 / (3.0 * fRef * fRef * fRef) + beta2 * log(fRef)) *
            etaInv +
            C1Int +
            C2Int * fRef,
            ifelse(
                fRef < fcutPar,
                (
                    -(alpha2 / fRef) +
                    (4.0 / 3.0) * (alpha3 * (fRef^(3.0 / 4.0))) +
                    alpha1 * fRef +
                    alpha4 * atan((fRef - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fRef,
                0.0,
            ),
        ),
    )

    phis = @. ifelse.(
        fgrid .< fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fgrid^(2.0 / 3.0) +
        PhiInspcoeffs.third * fgrid^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fgrid^(1.0 / 3.0) * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.log * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.min_third * fgrid^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fgrid^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fgrid +
        PhiInspcoeffs.min_four_thirds * fgrid^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fgrid^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fgrid +
            PhiInspcoeffs.four_thirds * fgrid^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fgrid^(5.0 / 3.0) +
            PhiInspcoeffs.two * fgrid * fgrid
        ) * etaInv,
        ifelse(
            fgrid < fMRDJoin,
            (beta1 * fgrid - beta3 / (3.0 * fgrid * fgrid * fgrid) + beta2 * log(fgrid)) * etaInv +
            C1Int +
            C2Int * fgrid,
            ifelse(
                fgrid < fcutPar,
                (
                    -(alpha2 / fgrid) +
                    (4.0 / 3.0) * (alpha3 * (fgrid^(3.0 / 4.0))) +
                    alpha1 * fgrid +
                    alpha4 * atan((fgrid - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fgrid,
                0.0,
            ),
        ),
    )
    if typeof(eta) == Float64
       phi = @. phis + ifelse(fgrid < fcutPar, -t0 * (fgrid - fRef) - phiRef, 0.0)
       if container !== nothing
        container .= [phi[i] for i in eachindex(phi)]
        end
    else 
        phi = @. phis + ifelse(fgrid .< fcutPar, -t0 * (fgrid - fRef) - phiRef, ForwardDiff.Dual{typeof(eta).parameters[1]}(0.,zeros(typeof(eta).parameters[3])...))
        if container !== nothing
            container .= [phi[i].value for i in eachindex(phi)]
        end
    end


    return phi
    #return @. phis + ifelse(fgrid < fcutPar, -t0 * (fgrid - fRef) - phiRef, 0.0)
end
""" helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Ampl(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing,
)
    # raise error if Lambda1 or Lambda2 are different from zero
    # if Lambda1 != 0.0 || Lambda2 != 0.0
    #     println("You requested Lambda1 and/or Lambda2 different from zero")
    #     println("Tidal parameters Lambda1 and Lambda2 are not supported in PhenomD")
    #     println("The code put them to zero")
    # end
    return Ampl(model, f, mc, eta, chi1, chi2, dL, fcutPar = fcutPar, fInsJoin_Ampl = fInsJoin_Ampl, GMsun_over_c3 = GMsun_over_c3, GMsun_over_c2_Gpc = GMsun_over_c2_Gpc, container = container)
end

"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomD(), f, mc, eta, chi1, chi2, dL)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Optional arguments:
-  `fInsJoin_Ampl`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The amplitude is dimensionless.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomD(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomD(), f, mc, eta, chi1, chi2, dL)
```
"""

function Ampl(model::PhenomD,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing,
    #container_jacobian = nothing,
)

    # """
    # Compute the amplitude of the GW as a function of frequency, given the events parameters.

    #  array f: Frequency grid on which the phase will be computed, in :math:`\\rm Hz`.
    #  dict(array, array, ...) kwargs: Dictionary with arrays containing the parameters of the events to compute the amplitude of, as in :py:data:`events`.
    #  GW amplitude for the chosen events evaluated on the frequency grid.
    # :rtype: array

    # """

    path=pwd()*"/useful_files/WFfiles/"
    QNMgrid_a = _readQNMgrid_a(path)
    QNMgrid_fring = _readQNMgrid_fring(path)
    QNMgrid_fdamp = _readQNMgrid_fdamp(path)
    # Useful quantities
    M = mc / (eta^(3.0 / 5.0))
    eta2 = eta * eta # This can speed up a bit, we call it multiple times
    pi2 = pi * pi
    chi12, chi22 = chi1 * chi1, chi2 * chi2

    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    # Compute ringdown and damping frequencies from interpolators
    fring = LinearInterpolation(QNMgrid_a, QNMgrid_fring)(aeff) / (1.0 - Erad)
    fdamp = LinearInterpolation(QNMgrid_a, QNMgrid_fdamp)(aeff) / (1.0 - Erad)
    # Compute coefficients gamma appearing in arXiv:1508.07253 eq. (19), the numerical coefficients are in Tab. 5
    xi2 = xi * xi
    gamma1 =
        0.006927402739328343 +
        0.03020474290328911 * eta +
        (
            0.006308024337706171 - 0.12074130661131138 * eta +
            0.26271598905781324 * eta2 +
            (
                0.0034151773647198794 - 0.10779338611188374 * eta +
                0.27098966966891747 * eta2
            ) * xi +
            (
                0.0007374185938559283 - 0.02749621038376281 * eta +
                0.0733150789135702 * eta2
            ) * xi2
        ) * xi
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi2
        ) * xi
    gamma3 =
        1.3081615607036106 - 0.005537729694807678 * eta +
        (
            -0.06782917938621007 - 0.6689834970767117 * eta +
            3.403147966134083 * eta2 +
            (-0.05296577374411866 - 0.9923793203111362 * eta + 4.820681208409587 * eta2) *
            xi +
            (
                -0.006134139870393713 - 0.38429253308696365 * eta +
                1.7561754421985984 * eta2
            ) * xi2
        ) * xi
    # Compute fpeak, from arXiv:1508.07253 eq. (20), we remove the square root term in case it is complex
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs(fring - (fdamp * gamma3) / gamma2),
        fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2 * gamma2)) * gamma3) / gamma2,
    )
    # Compute coefficients rho appearing in arXiv:1508.07253 eq. (30), the numerical coefficients are in Tab. 5
    rho1 =
        3931.8979897196696 - 17395.758706812805 * eta +
        (
            3132.375545898835 + 343965.86092361377 * eta - 1.2162565819981997e6 * eta2 +
            (-70698.00600428853 + 1.383907177859705e6 * eta - 3.9662761890979446e6 * eta2) *
            xi +
            (-60017.52423652596 + 803515.1181825735 * eta - 2.091710365941658e6 * eta2) *
            xi2
        ) * xi
    rho2 =
        -40105.47653771657 +
        112253.0169706701 * eta +
        (
            23561.696065836168 - 3.476180699403351e6 * eta +
            1.137593670849482e7 * eta2 +
            (754313.1127166454 - 1.308476044625268e7 * eta + 3.6444584853928134e7 * eta2) *
            xi +
            (596226.612472288 - 7.4277901143564405e6 * eta + 1.8928977514040343e7 * eta2) *
            xi2
        ) * xi
    rho3 =
        83208.35471266537 - 191237.7264145924 * eta +
        (
            -210916.2454782992 + 8.71797508352568e6 * eta - 2.6914942420669552e7 * eta2 +
            (
                -1.9889806527362722e6 + 3.0888029960154563e7 * eta -
                8.390870279256162e7 * eta2
            ) * xi +
            (
                -1.4535031953446497e6 + 1.7063528990822166e7 * eta -
                4.2748659731120914e7 * eta2
            ) * xi2
        ) * xi
    # Compute coefficients delta appearing in arXiv:1508.07253 eq. (21)
    f1Interm = fInsJoin_Ampl
    f3Interm = fpeak
    dfInterm = 0.5 * (f3Interm - f1Interm)
    f2Interm = f1Interm + dfInterm
    # First write the inspiral coefficients, we put them in a dictionary and label with the power in front of which they appear
    amp0 = sqrt(2.0 * eta / 3.0) * (pi^(-1.0 / 6.0))
    Acoeffs = AcoeffsStructure(
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / 672.0,
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0,
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) * (pi^(4.0 / 3.0))
        ) / 8.128512e6,
        (
            (pi^(5.0 / 3.0)) * (
                chi2 * (
                    -285197.0 * (-1.0 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1.0 + 4.0 * eta) * pi
            )
        ) / 32256.0,
        -(
            (pi^2.0) * (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368.0 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320.0 * pi2
                )
            )
        ) / 6.0085960704e10,
        rho1,
        rho2,
        rho3,
    )

    # v1 is the inspiral model evaluated at f1Interm
    v1 =
        1.0 +
        (f1Interm^(2.0 / 3.0)) * Acoeffs.two_thirds +
        (f1Interm^(4.0 / 3.0)) * Acoeffs.four_thirds +
        (f1Interm^(5.0 / 3.0)) * Acoeffs.five_thirds +
        (f1Interm^(7.0 / 3.0)) * Acoeffs.seven_thirds +
        (f1Interm^(8.0 / 3.0)) * Acoeffs.eight_thirds +
        f1Interm *
        (Acoeffs.one + f1Interm * Acoeffs.two + f1Interm * f1Interm * Acoeffs.three)
    # d1 is the derivative of the inspiral model evaluated at f1
    d1 =
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / (1008.0 * (f1Interm^(1.0 / 3.0))) +
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0 +
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) *
            (f1Interm^(1.0 / 3.0)) *
            (pi^(4.0 / 3.0))
        ) / 6.096384e6 +
        (
            5.0 *
            (f1Interm^(2.0 / 3.0)) *
            (pi^(5.0 / 3.0)) *
            (
                chi2 * (
                    -285197.0 * (-1 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1 + 4 * eta) * pi
            )
        ) / 96768.0 -
        (
            f1Interm *
            pi2 *
            (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320 * pi2
                )
            )
        ) / 3.0042980352e10 +
        (7.0 / 3.0) * (f1Interm^(4.0 / 3.0)) * rho1 +
        (8.0 / 3.0) * (f1Interm^(5.0 / 3.0)) * rho2 +
        3.0 * (f1Interm * f1Interm) * rho3
    # v3 is the merger-ringdown model (eq. (19) of arXiv:1508.07253) evaluated at f3
    v3 =
        exp(-(f3Interm - fring) * gamma2 / (fdamp * gamma3)) * (fdamp * gamma3 * gamma1) /
        ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
    # d2 is the derivative of the merger-ringdown model evaluated at f3
    d2 =
        (
            (-2.0 * fdamp * (f3Interm - fring) * gamma3 * gamma1) /
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3) -
            (gamma2 * gamma1)
        ) / (
            exp((f3Interm - fring) * gamma2 / (fdamp * gamma3)) *
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
        )
    # v2 is the value of the amplitude evaluated at f2. They come from the fit of the collocation points in the intermediate region
    v2 =
        0.8149838730507785 +
        2.5747553517454658 * eta +
        (
            1.1610198035496786 - 2.3627771785551537 * eta +
            6.771038707057573 * eta2 +
            (0.7570782938606834 - 2.7256896890432474 * eta + 7.1140380397149965 * eta2) *
            xi +
            (0.1766934149293479 - 0.7978690983168183 * eta + 2.1162391502005153 * eta2) *
            xi2
        ) * xi
    # Now some definitions to speed up
    f1 = f1Interm
    f2 = f2Interm
    f3 = f3Interm
    f12 = f1Interm * f1Interm
    f13 = f1Interm * f12
    f14 = f1Interm * f13
    f15 = f1Interm * f14
    f22 = f2Interm * f2Interm
    f23 = f2Interm * f22
    f24 = f2Interm * f23
    f32 = f3Interm * f3Interm
    f33 = f3Interm * f32
    f34 = f3Interm * f33
    f35 = f3Interm * f34
    # Finally conpute the deltas
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2.0 * d1 * f13 * f22 * f33 - 2.0 * d2 * f13 * f22 * f33 -
            d1 * f12 * f23 * f33 + d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 -
            d1 * f13 * f2 * f34 - d1 * f12 * f22 * f34 +
            2.0 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4.0 * f12 * f23 * f32 * v1 -
            3.0 * f1 * f24 * f32 * v1 - 8.0 * f12 * f22 * f33 * v1 +
            4.0 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4.0 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2.0 * f23 * f34 * v1 - 2.0 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3.0 * f14 * f33 * v2 -
            3.0 * f13 * f34 * v2 + f12 * f35 * v2 - f15 * f22 * v3 +
            2.0 * f14 * f23 * v3 - f13 * f24 * v3 + 2.0 * f15 * f2 * f3 * v3 -
            f14 * f22 * f3 * v3 - 4.0 * f13 * f23 * f3 * v3 +
            3.0 * f12 * f24 * f3 * v3 - 4.0 * f14 * f2 * f32 * v3 +
            8.0 * f13 * f22 * f32 * v3 - 4.0 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2 * d1 * f13 * f22 * f33 - 2 * d2 * f13 * f22 * f33 - d1 * f12 * f23 * f33 +
            d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 - d1 * f13 * f2 * f34 -
            d1 * f12 * f22 * f34 +
            2 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4 * f12 * f23 * f32 * v1 -
            3 * f1 * f24 * f32 * v1 - 8 * f12 * f22 * f33 * v1 +
            4 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2 * f23 * f34 * v1 - 2 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3 * f14 * f33 * v2 - 3 * f13 * f34 * v2 +
            f12 * f35 * v2 - f15 * f22 * v3 + 2 * f14 * f23 * v3 - f13 * f24 * v3 +
            2 * f15 * f2 * f3 * v3 - f14 * f22 * f3 * v3 - 4 * f13 * f23 * f3 * v3 +
            3 * f12 * f24 * f3 * v3 - 4 * f14 * f2 * f32 * v3 +
            8 * f13 * f22 * f32 * v3 - 4 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta1 = -(
        (
            -(d2 * f15 * f22) + 2.0 * d2 * f14 * f23 - d2 * f13 * f24 -
            d2 * f14 * f22 * f3 +
            2.0 * d1 * f13 * f23 * f3 +
            2.0 * d2 * f13 * f23 * f3 - 2 * d1 * f12 * f24 * f3 - d2 * f12 * f24 * f3 +
            d2 * f15 * f32 - 3 * d1 * f13 * f22 * f32 - d2 * f13 * f22 * f32 +
            2 * d1 * f12 * f23 * f32 - 2 * d2 * f12 * f23 * f32 +
            d1 * f1 * f24 * f32 +
            2 * d2 * f1 * f24 * f32 - d2 * f14 * f33 +
            d1 * f12 * f22 * f33 +
            3 * d2 * f12 * f22 * f33 - 2 * d1 * f1 * f23 * f33 -
            2 * d2 * f1 * f23 * f33 +
            d1 * f24 * f33 +
            d1 * f13 * f34 +
            d1 * f1 * f22 * f34 - 2 * d1 * f23 * f34 - d1 * f12 * f35 + d1 * f22 * f35 -
            8 * f12 * f23 * f3 * v1 +
            6 * f1 * f24 * f3 * v1 +
            12 * f12 * f22 * f32 * v1 - 8 * f1 * f23 * f32 * v1 - 4 * f12 * f34 * v1 +
            2 * f1 * f35 * v1 +
            2 * f15 * f3 * v2 - 4 * f14 * f32 * v2 + 4 * f12 * f34 * v2 -
            2 * f1 * f35 * v2 - 2 * f15 * f3 * v3 + 8 * f12 * f23 * f3 * v3 -
            6 * f1 * f24 * f3 * v3 + 4 * f14 * f32 * v3 - 12 * f12 * f22 * f32 * v3 +
            8 * f1 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta2 = -(
        (
            d2 * f15 * f2 - d1 * f13 * f23 - 3 * d2 * f13 * f23 +
            d1 * f12 * f24 +
            2.0 * d2 * f12 * f24 - d2 * f15 * f3 + d2 * f14 * f2 * f3 -
            d1 * f12 * f23 * f3 +
            d2 * f12 * f23 * f3 +
            d1 * f1 * f24 * f3 - d2 * f1 * f24 * f3 - d2 * f14 * f32 +
            3 * d1 * f13 * f2 * f32 +
            d2 * f13 * f2 * f32 - d1 * f1 * f23 * f32 + d2 * f1 * f23 * f32 -
            2 * d1 * f24 * f32 - d2 * f24 * f32 - 2 * d1 * f13 * f33 +
            2 * d2 * f13 * f33 - d1 * f12 * f2 * f33 - 3 * d2 * f12 * f2 * f33 +
            3 * d1 * f23 * f33 +
            d2 * f23 * f33 +
            d1 * f12 * f34 - d1 * f1 * f2 * f34 + d1 * f1 * f35 - d1 * f2 * f35 +
            4 * f12 * f23 * v1 - 3 * f1 * f24 * v1 + 4 * f1 * f23 * f3 * v1 -
            3 * f24 * f3 * v1 - 12 * f12 * f2 * f32 * v1 +
            4 * f23 * f32 * v1 +
            8 * f12 * f33 * v1 - f1 * f34 * v1 - f35 * v1 - f15 * v2 - f14 * f3 * v2 +
            8 * f13 * f32 * v2 - 8 * f12 * f33 * v2 +
            f1 * f34 * v2 +
            f35 * v2 +
            f15 * v3 - 4 * f12 * f23 * v3 +
            3 * f1 * f24 * v3 +
            f14 * f3 * v3 - 4 * f1 * f23 * f3 * v3 + 3 * f24 * f3 * v3 -
            8 * f13 * f32 * v3 + 12 * f12 * f2 * f32 * v3 - 4 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta3 = -(
        (
            -2.0 * d2 * f14 * f2 + d1 * f13 * f22 + 3 * d2 * f13 * f22 - d1 * f1 * f24 -
            d2 * f1 * f24 + 2 * d2 * f14 * f3 - 2.0 * d1 * f13 * f2 * f3 -
            2 * d2 * f13 * f2 * f3 + d1 * f12 * f22 * f3 - d2 * f12 * f22 * f3 +
            d1 * f24 * f3 +
            d2 * f24 * f3 +
            d1 * f13 * f32 - d2 * f13 * f32 - 2 * d1 * f12 * f2 * f32 +
            2 * d2 * f12 * f2 * f32 +
            d1 * f1 * f22 * f32 - d2 * f1 * f22 * f32 + d1 * f12 * f33 -
            d2 * f12 * f33 +
            2 * d1 * f1 * f2 * f33 +
            2 * d2 * f1 * f2 * f33 - 3 * d1 * f22 * f33 - d2 * f22 * f33 -
            2 * d1 * f1 * f34 + 2 * d1 * f2 * f34 - 4 * f12 * f22 * v1 +
            2 * f24 * v1 +
            8 * f12 * f2 * f3 * v1 - 4 * f1 * f22 * f3 * v1 - 4 * f12 * f32 * v1 +
            8 * f1 * f2 * f32 * v1 - 4 * f22 * f32 * v1 - 4 * f1 * f33 * v1 +
            2 * f34 * v1 +
            2 * f14 * v2 - 4 * f13 * f3 * v2 + 4 * f1 * f33 * v2 - 2 * f34 * v2 -
            2 * f14 * v3 + 4 * f12 * f22 * v3 - 2 * f24 * v3 + 4 * f13 * f3 * v3 -
            8 * f12 * f2 * f3 * v3 +
            4 * f1 * f22 * f3 * v3 +
            4 * f12 * f32 * v3 - 8 * f1 * f2 * f32 * v3 + 4 * f22 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta4 = -(
        (
            d2 * f13 * f2 - d1 * f12 * f22 - 2 * d2 * f12 * f22 +
            d1 * f1 * f23 +
            d2 * f1 * f23 - d2 * f13 * f3 +
            2.0 * d1 * f12 * f2 * f3 +
            d2 * f12 * f2 * f3 - d1 * f1 * f22 * f3 + d2 * f1 * f22 * f3 -
            d1 * f23 * f3 - d2 * f23 * f3 - d1 * f12 * f32 + d2 * f12 * f32 -
            d1 * f1 * f2 * f32 - 2 * d2 * f1 * f2 * f32 +
            2 * d1 * f22 * f32 +
            d2 * f22 * f32 +
            d1 * f1 * f33 - d1 * f2 * f33 + 3 * f1 * f22 * v1 - 2 * f23 * v1 -
            6 * f1 * f2 * f3 * v1 +
            3 * f22 * f3 * v1 +
            3 * f1 * f32 * v1 - f33 * v1 - f13 * v2 + 3 * f12 * f3 * v2 -
            3 * f1 * f32 * v2 +
            f33 * v2 +
            f13 * v3 - 3 * f1 * f22 * v3 + 2 * f23 * v3 - 3 * f12 * f3 * v3 +
            6 * f1 * f2 * f3 * v3 - 3 * f22 * f3 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )

    # Defined as in LALSimulation - LALSimIMRPhenomD.c line 332. Final units are correctly Hz^-1
    Overallamp =
        2.0 * sqrt(5.0 / (64.0 * pi)) * M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL

    amplitudeIMR = @. ifelse(
        fgrid < fInsJoin_Ampl,
        1.0 +
        (fgrid^(2.0 / 3.0)) * Acoeffs.two_thirds +
        (fgrid^(4.0 / 3.0)) * Acoeffs.four_thirds +
        (fgrid^(5.0 / 3.0)) * Acoeffs.five_thirds +
        (fgrid^(7.0 / 3.0)) * Acoeffs.seven_thirds +
        (fgrid^(8.0 / 3.0)) * Acoeffs.eight_thirds +
        fgrid * (Acoeffs.one + fgrid * Acoeffs.two + fgrid * fgrid * Acoeffs.three),
        ifelse(
            fgrid < fpeak,
            delta0 +
            fgrid * delta1 +
            fgrid * fgrid * (delta2 + fgrid * delta3 + fgrid * fgrid * delta4),
            ifelse(
                fgrid < fcutPar,
                exp(-(fgrid - fring) * gamma2 / (fdamp * gamma3)) *
                (fdamp * gamma3 * gamma1) /
                ((fgrid - fring) * (fgrid - fring) + fdamp * gamma3 * fdamp * gamma3),
                0.0,
            ),
        ),
    )
    if typeof(eta) == Float64
        ampl = Overallamp * amp0 .* (fgrid .^ (-7.0 / 6.0)) .* amplitudeIMR
        if container !== nothing
            container .= [ampl[i] for i in eachindex(ampl)]
         end
    else 
         #phi = @. phis + ifelse(fgrid .< fcutPar, -t0 * (fgrid - fRef) - phiRef, ForwardDiff.Dual{typeof(eta).parameters[1]}(0.,zeros(typeof(eta).parameters[3])...))
        ampl = Overallamp * amp0 .* (fgrid .^ (-7.0 / 6.0)) .* amplitudeIMR

        if container !== nothing
            container .= [ampl[i].value for i in eachindex(ampl)]
        end
        # if container_jacobian !== nothing
        #     container_jacobian .= [ampl[i] for i in eachindex(ampl)]
        # end
    end
    # ampl = Overallamp * amp0 .* (fgrid .^ (-7.0 / 6.0)) .* amplitudeIMR

    # if container !== nothing
    #     container .= [ampl[i].value for i in eachindex(ampl)]
    # end

    return ampl
end

"""
Compute the spin of the final object, as in LALSimIMRPhenomD_internals.c line 161 and 142, which is taken from `arXiv:1508.07250 <https://arxiv.org/abs/1508.07250>`_ eq. (3.6).
Valid for all the models considered in this code.


#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.

#### Return:
-  (float) spin of the final object.

"""


function _finalspin(model::Model, eta, chi1, chi2)

    eta2 = eta * eta
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    s = (m1 * m1 * chi1 + m2 * m2 * chi2)

    af1 =
        eta * (
            3.4641016151377544 - 4.399247300629289 * eta + 9.397292189321194 * eta2 -
            13.180949901606242 * eta2 * eta
        )
    af2 =
        eta * (
            s * (
                (1.0 / eta - 0.0850917821418767 - 5.837029316602263 * eta) +
                (0.1014665242971878 - 2.0967746996832157 * eta) * s
            )
        )
    af3 =
        eta * (
            s * (
                (-1.3546806617824356 + 4.108962025369336 * eta) * s * s +
                (-0.8676969352555539 + 2.064046835273906 * eta) * s * s * s
            )
        )
    return af1 + af2 + af3
end

"""
Compute the total radiated energy, as in `arXiv:1508.07250 <https://arxiv.org/abs/1508.07250>`_ eq. (3.7) and (3.8).
Valid for all the models considered in this code.

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.

#### Return:
-  (float) total radiated energy.

"""

function _radiatednrg(model::Model, eta, chi1, chi2)


    eta2 = eta * eta
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    m1 = 0.5 * (1.0 + Seta)
    m2 = 0.5 * (1.0 - Seta)
    s = (m1 * m1 * chi1 + m2 * m2 * chi2) / (m1 * m1 + m2 * m2)

    EradNS =
        eta * (
            0.055974469826360077 + 0.5809510763115132 * eta - 0.9606726679372312 * eta2 +
            3.352411249771192 * eta2 * eta
        )

    return (
        EradNS * (
            1.0 +
            (
                -0.0030302335878845507 - 2.0066110851351073 * eta +
                7.7050567802399215 * eta2
            ) * s
        )
    ) / (
        1.0 +
        (-0.6714403054720589 - 1.4756929437702908 * eta + 7.304676214885011 * eta2) * s
    )
end


"""
Compute the time to coalescence (in seconds) as a function of frequency in `Hz`. Used when including the Earth motion effect.

We use the expression in `arXiv:0907.0700 <https://arxiv.org/abs/0907.0700>`_ eq. (3.8b).

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency grid for a single event on which the time to coalescence will be computed, in `Hz`
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.

#### Return:
-  (array) time to coalescence for the chosen event evaluated on the frequency grid, in seconds.

"""

function _tau_star(model::Model, f, mc, eta; GMsun_over_c3 = uc.GMsun_over_c3)


    Mtot_sec = mc * GMsun_over_c3 / (eta^(0.6))
    v = (pi * Mtot_sec .* f) .^ (1.0 / 3.0)
    eta2 = eta * eta

    OverallFac = @. 5.0 / 256 * Mtot_sec / (eta * (v .^ 8.0))
    v2 = v .* v
    v3 = v2 .* v
    v4 = v3 .* v
    v5 = v4 .* v
    v6 = v5 .* v

    t05 =
        1.0 .+ (743.0 / 252.0 + 11.0 / 3.0 * eta) .* v2 .- 32.0 / 5.0 * pi .* v3 .+
        (3058673.0 / 508032.0 + 5429.0 / 504.0 * eta + 617.0 / 72.0 * eta2) .* v4 .-
        (7729.0 / 252.0 - 13.0 / 3.0 * eta) * pi .* v5
    t6 = 
        @. (
            -10052469856691.0 / 23471078400.0 +
            128.0 / 3.0 * pi^2 +
            6848.0 / 105.0 * MathConstants.eulergamma +
            (3147553127.0 / 3048192.0 - 451.0 / 12.0 * pi^2) * eta -
            15211.0 / 1728.0 * eta2 +
            25565.0 / 1296.0 * eta2 * eta +
            3424.0 / 105.0 * log(16.0 .* v2)
        ) * (v6)
    t7 = 
        @. (-15419335.0 / 127008.0 - 75703.0 / 756.0 * eta + 14809.0 / 378.0 * eta2) * pi .*
        v6 .* v

    return OverallFac .* (t05 .+ t6 .+ t7)
end

function _fcut(model::Model, mc, eta, Lambda1, Lambda2; fcutPar = 0.2, GMsun_over_c3 = uc.GMsun_over_c3)
    return _fcut(model, mc, eta, fcutPar = fcutPar, GMsun_over_c3 = GMsun_over_c3)
end
"""
Compute the cut frequency of the waveform as a function of the events parameters, in `Hz`.
Valid for TaylorF2, PhenomD, PhenomHM. Not for PhenomD_NRTidal.
We use the expression in `arXiv:0907.0700 <https://arxiv.org/abs/0907.0700>`_ eq. (3.8b).

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.

#### Optional arguments:
-  `fcutPar`: The cut frequency factor of the waveform as an adimensional frequency (Mf). Default is 0.2.

#### Return:
-  (float) Cut frequency of the waveform, in Hz.

"""
function _fcut(model::Model, mc, eta; fcutPar = 0.2, GMsun_over_c3 = uc.GMsun_over_c3)

    if typeof(model)== PhenomD_NRTidal
        error("You need to provide also Lambda1 and Lambda2 since fcut depends on them")
    end
    if typeof(model) == TaylorF2
        return 1. / ( 6 * pi * sqrt(6.) * GMsun_over_c3 )
    end

    return fcutPar / (mc * GMsun_over_c3 / (eta^(0.6)))

end


##############################################################################
# IMRPhenomD_NRTidalv2 WAVEFORM
##############################################################################
"""
IMRPhenomD_NRTidal waveform model.

Relevant references:
    [1] `arXiv:1508.07250 <https://arxiv.org/abs/1508.07250>`_
    
    [2] `arXiv:1508.07253 <https://arxiv.org/abs/1508.07253>`_
    
    [3] `arXiv:1905.06011 <https://arxiv.org/abs/1905.06011>`_

 kwargs: Optional arguments to be passed to the parent class :py:class:`WaveFormModel`, such as ``is_chi1chi2``.
    
"""
# All is taken from LALSimulation and arXiv:1508.07250, arXiv:1508.07253, arXiv:1905.06011

"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomD_NRTidal(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `Lambda1`: Dimensionless tidal deformability of the first NS.
-  `Lambda2`: Dimensionless tidal deformability of the second NS.
#### Optional arguments:
-  `fInsJoin_PHI`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The phase is given in radians.

#### Example:
```julia
    mc = 3.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    Lambda1 = 1000.
    Lambda2 = 2000.
    fcut = _fcut(PhenomD_NRTidal(), mc, eta, Lambda1, Lambda2)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Phi(PhenomD_NRTidal(), f, mc, eta, dL, chi1, chi2, Lambda1, Lambda2)
```
"""

function Phi(model::PhenomD_NRTidal,
    f,
    mc,
    eta,
    chi1,
    chi2,
    Lambda1,
    Lambda2;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)
    # # print that if Lambda1 and or Lambda2 are zero one should resolve to other waveforms
    # if Lambda1 == 0. && Lambda2 ==0.
    #     println("Lambda1 and Lambda2 are zero, please use other waveforms, valid for BBH")
    # elseif Lambda1 == 0. || Lambda2 == 0.
    #     println("Lambda1 or Lambda2 is zero, please use PhenomNSBH")
    # end

    path=pwd()*"/useful_files/WFfiles/"
    QNMgrid_a = _readQNMgrid_a(path)
    QNMgrid_fring = _readQNMgrid_fring(path)
    QNMgrid_fdamp = _readQNMgrid_fdamp(path)

    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta

    pi2 = pi * pi

    fInsJoin = fInsJoin_PHI

        # A non-zero tidal deformability induces a quadrupole moment (for BBH it is 1).
        # The relation between the two is given in arxiv:1608.02582 eq. (15) with coefficients from third row of Table I
        # We also extend the range to 0 <= Lam < 1, as done in LALSimulation in LALSimUniversalRelations.c line 123
    QuadMon1 = ifelse(Lambda1 < 1., 1. + Lambda1*(0.427688866723244 + Lambda1*(-0.324336526985068 + Lambda1*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda1) + 0.04812 * log(Lambda1)^2 -4.283e-3 * log(Lambda1)^3 + 1.245e-4 * log(Lambda1)^4))
    QuadMon2 = ifelse(Lambda2 < 1., 1. + Lambda2*(0.427688866723244 + Lambda2*(-0.324336526985068 + Lambda2*0.1107439432180572)), exp(0.1940 + 0.09163 * log(Lambda2) + 0.04812 * log(Lambda2)^2 -4.283e-3 * log(Lambda2)^3 + 1.245e-4 * log(Lambda2)^4))
    

    chi12, chi22 = chi1*chi1, chi2*chi2
    chi1dotchi2  = chi1*chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)       
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    chi_s2, chi_a2 = chi_s*chi_s, chi_a*chi_a
    chi_sdotchi_a  = chi_s*chi_a
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M*GMsun_over_c3.*f
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = - 1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    # Compute ringdown and damping frequencies from interpolators
    fring = LinearInterpolation(QNMgrid_a, QNMgrid_fring)(aeff) / (1.0 - Erad)
    fdamp = LinearInterpolation(QNMgrid_a, QNMgrid_fdamp)(aeff) / (1.0 - Erad)

    # Compute sigma coefficients appearing in arXiv:1508.07253 eq. (28)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    sigma1 =
        2096.551999295543 +
        1463.7493168261553 * eta +
        (
            1312.5493286098522 + 18307.330017082117 * eta - 43534.1440746107 * eta2 +
            (-833.2889543511114 + 32047.31997183187 * eta - 108609.45037520859 * eta2) *
            xi +
            (452.25136398112204 + 8353.439546391714 * eta - 44531.3250037322 * eta2) * xi2
        ) * xi
    sigma2 =
        -10114.056472621156 - 44631.01109458185 * eta +
        (
            -6541.308761668722 - 266959.23419307504 * eta +
            686328.3229317984 * eta2 +
            (3405.6372187679685 - 437507.7208209015 * eta + 1.6318171307344697e6 * eta2) *
            xi +
            (-7462.648563007646 - 114585.25177153319 * eta + 674402.4689098676 * eta2) * xi2
        ) * xi
    sigma3 =
        22933.658273436497 +
        230960.00814979506 * eta +
        (
            14961.083974183695 + 1.1940181342318142e6 * eta - 3.1042239693052764e6 * eta2 +
            (-3038.166617199259 + 1.8720322849093592e6 * eta - 7.309145012085539e6 * eta2) *
            xi +
            (42738.22871475411 + 467502.018616601 * eta - 3.064853498512499e6 * eta2) * xi2
        ) * xi
    sigma4 =
        -14621.71522218357 - 377812.8579387104 * eta +
        (
            -9608.682631509726 - 1.7108925257214056e6 * eta +
            4.332924601416521e6 * eta2 +
            (
                -22366.683262266528 - 2.5019716386377467e6 * eta +
                1.0274495902259542e7 * eta2
            ) * xi +
            (-85360.30079034246 - 570025.3441737515 * eta + 4.396844346849777e6 * eta2) *
            xi2
        ) * xi

    # Compute beta coefficients appearing in arXiv:1508.07253 eq. (16)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    beta1 =
        97.89747327985583 - 42.659730877489224 * eta +
        (
            153.48421037904913 - 1417.0620760768954 * eta +
            2752.8614143665027 * eta2 +
            (138.7406469558649 - 1433.6585075135881 * eta + 2857.7418952430758 * eta2) *
            xi +
            (41.025109467376126 - 423.680737974639 * eta + 850.3594335657173 * eta2) * xi2
        ) * xi
    beta2 =
        -3.282701958759534 - 9.051384468245866 * eta +
        (
            -12.415449742258042 + 55.4716447709787 * eta - 106.05109938966335 * eta2 +
            (-11.953044553690658 + 76.80704618365418 * eta - 155.33172948098394 * eta2) *
            xi +
            (-3.4129261592393263 + 25.572377569952536 * eta - 54.408036707740465 * eta2) *
            xi2
        ) * xi
    beta3 =
        -0.000025156429818799565 +
        0.000019750256942201327 * eta +
        (
            -0.000018370671469295915 +
            0.000021886317041311973 * eta +
            0.00008250240316860033 * eta2 +
            (
                7.157371250566708e-6 - 0.000055780000112270685 * eta +
                0.00019142082884072178 * eta2
            ) * xi +
            (
                5.447166261464217e-6 - 0.00003220610095021982 * eta +
                0.00007974016714984341 * eta2
            ) * xi2
        ) * xi

    # Compute alpha coefficients appearing in arXiv:1508.07253 eq. (14)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    alpha1 =
        43.31514709695348 +
        638.6332679188081 * eta +
        (
            -32.85768747216059 + 2415.8938269370315 * eta - 5766.875169379177 * eta2 +
            (-61.85459307173841 + 2953.967762459948 * eta - 8986.29057591497 * eta2) * xi +
            (-21.571435779762044 + 981.2158224673428 * eta - 3239.5664895930286 * eta2) *
            xi2
        ) * xi
    alpha2 =
        -0.07020209449091723 - 0.16269798450687084 * eta +
        (
            -0.1872514685185499 + 1.138313650449945 * eta - 2.8334196304430046 * eta2 +
            (-0.17137955686840617 + 1.7197549338119527 * eta - 4.539717148261272 * eta2) *
            xi +
            (-0.049983437357548705 + 0.6062072055948309 * eta - 1.682769616644546 * eta2) *
            xi2
        ) * xi
    alpha3 =
        9.5988072383479 - 397.05438595557433 * eta +
        (
            16.202126189517813 - 1574.8286986717037 * eta +
            3600.3410843831093 * eta2 +
            (27.092429659075467 - 1786.482357315139 * eta + 5152.919378666511 * eta2) * xi +
            (11.175710130033895 - 577.7999423177481 * eta + 1808.730762932043 * eta2) * xi2
        ) * xi
    alpha4 =
        -0.02989487384493607 +
        1.4022106448583738 * eta +
        (
            -0.07356049468633846 +
            0.8337006542278661 * eta +
            0.2240008282397391 * eta2 +
            (-0.055202870001177226 + 0.5667186343606578 * eta + 0.7186931973380503 * eta2) *
            xi +
            (
                -0.015507437354325743 +
                0.15750322779277187 * eta +
                0.21076815715176228 * eta2
            ) * xi2
        ) * xi
    alpha5 =
        0.9974408278363099 - 0.007884449714907203 * eta +
        (
            -0.059046901195591035 + 1.3958712396764088 * eta - 4.516631601676276 * eta2 +
            (-0.05585343136869692 + 1.7516580039343603 * eta - 5.990208965347804 * eta2) *
            xi +
            (-0.017945336522161195 + 0.5965097794825992 * eta - 2.0608879367971804 * eta2) *
            xi2
        ) * xi

    # Compute the TF2 phase coefficients and put them in a dictionary (spin effects are included up to 3.5PN)
    TF2OverallAmpl = 3 / (128.0 * eta)

    TF2_5coeff_tmp =
        38645.0 * pi / 756.0 - 65.0 * pi * eta / 9.0 - (
            (732985.0 / 2268.0 - 24260.0 * eta / 81.0 - 340.0 * eta2 / 9.0) * chi_s +
            (732985.0 / 2268.0 + 140.0 * eta / 9.0) * Seta * chi_a
        ) #variable to be used later
    TF2_6coeff_tmp =
        11583.231236531 / 4.694215680 - 640.0 / 3.0 * pi2 -
        684.8 / 2.1 * MathConstants.eulergamma +
        eta * (-15737.765635 / 3.048192 + 225.5 / 1.2 * pi2) +
        eta2 * 76.055 / 1.728 - eta2 * eta * 127.825 / 1.296 - log(4.0) * 684.8 / 2.1 +
        pi * chi1 * m1ByM * (1490.0 / 3.0 + m1ByM * 260.0) +
        pi * chi2 * m2ByM * (1490.0 / 3.0 + m2ByM * 260.0) +
        (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) *
        m1ByM^2 *
        QuadMon1 *
        chi12 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 ) *
        m1ByM^2 *
        chi12 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) *
        m2ByM^2 *
        QuadMon2 *
        chi22 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 ) *
        m2ByM^2 *
        chi22

    TF2coeffs = TF2coeffsStructure(
        1.0,
        0.0,
        3715.0 / 756.0 + (55.0 * eta) / 9.0,
        -16.0 * pi +
        (113.0 * Seta * chi_a) / 3.0 +
        (113.0 / 3.0 - (76.0 * eta) / 3.0) * chi_s,
        5.0 * (3058.673 / 7.056 + 5429.0 / 7.0 * eta + 617.0 * eta2) / 72.0 +
        247.0 / 4.8 * eta * chi1dotchi2 - 721.0 / 4.8 * eta * chi1dotchi2 +
        (-720.0 / 9.6 * QuadMon1 + 1.0 / 9.6) * m1ByM^2  * chi12 +
        (-720.0 / 9.6 * QuadMon2 + 1.0 / 9.6) * m2ByM^2  * chi22 +
        (240.0 / 9.6 * QuadMon1 - 7.0 / 9.6) * m1ByM^2  * chi12 +
        (240.0 / 9.6 * QuadMon2 - 7.0 / 9.6) * m2ByM^2  * chi22,
        TF2_5coeff_tmp,
        TF2_5coeff_tmp * 3.0,
        TF2_6coeff_tmp - (
            (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 )
            ) *
            m1ByM *
            m1ByM *
            chi12 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 )
            ) *
            m2ByM *
            m2ByM *
            chi22
        ),
        -6848.0 / 21.0,
        77096675.0 * pi / 254016.0 + 378515.0 * pi * eta / 1512.0 -
        74045.0 * pi * eta2 / 756.0 +
        (
            -25150083775.0 / 3048192.0 + 10566655595.0 * eta / 762048.0 -
            1042165.0 * eta2 / 3024.0 + 5345.0 * eta2 * eta / 36.0
        ) * chi_s +
        Seta * (
            (
                -25150083775.0 / 3048192.0 + 26804935.0 * eta / 6048.0 -
                1985.0 * eta2 / 48.0
            ) * chi_a
        ),
    )


    PhiInspcoeffs = PhiInspcoeffsStructure(
        TF2coeffs.five * TF2OverallAmpl,
        TF2coeffs.seven * TF2OverallAmpl * (pi^(2.0 / 3.0)),
        TF2coeffs.six * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.six_log * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.five_log * TF2OverallAmpl,
        TF2coeffs.four * TF2OverallAmpl * (pi^(-1.0 / 3.0)),
        TF2coeffs.three * TF2OverallAmpl * (pi^(-2.0 / 3.0)),
        TF2coeffs.two * TF2OverallAmpl / pi,
        TF2coeffs.one * TF2OverallAmpl * (pi^(-4.0 / 3.0)),
        TF2coeffs.zero * TF2OverallAmpl * (pi^(-5.0 / 3.0)),
        sigma1,
        sigma2 * 0.75,
        sigma3 * 0.6,
        sigma4 * 0.5,
    )

    #initial_phasing 
    #log
    #Now compute the coefficients to align the three parts

    fMRDJoin = 0.5 * fring

    # First the Inspiral - Intermediate: we compute C1Int and C2Int coeffs
    # Equations to solve for to get C(1) continuous join
    # PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    # Joining at fInsJoin
    # PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    # PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    # This is the first derivative wrt f of the inspiral phase computed at fInsJoin, first add the PN contribution and then the higher order calibrated terms

    #DPhiIns = (2.0*TF2coeffs['seven']*TF2OverallAmpl*((pi*fInsJoin)^(7. /3.)) + (TF2coeffs['six']*TF2OverallAmpl + TF2coeffs['six_log']*TF2OverallAmpl * (1.0 + log(pi*fInsJoin)/3.))*((pi*fInsJoin)^(2.)) + TF2coeffs['five_log']*TF2OverallAmpl*((pi*fInsJoin)^(5. /3.)) - TF2coeffs['four']*TF2OverallAmpl*((pi*fInsJoin)^(4. /3.)) - 2. *TF2coeffs['three']*TF2OverallAmpl*(pi*fInsJoin) - 3. *TF2coeffs['two']*TF2OverallAmpl*((pi*fInsJoin)^(2. /3.)) - 4. *TF2coeffs['one']*TF2OverallAmpl*((pi*fInsJoin)^(1. /3.)) - 5. *TF2coeffs['zero']*TF2OverallAmpl)*pi/(3. *((pi*fInsJoin)^(8. /3.)))
    #DPhiIns = (2.0*TF2coeffs['seven']*TF2OverallAmpl*((pi*fInsJoin)^(7. /3.)) + (TF2coeffs['six']*TF2OverallAmpl + TF2coeffs['six_log']*TF2OverallAmpl * (1.0 + log(pi*fInsJoin)/3.))*((pi*fInsJoin)^(2.)) + TF2coeffs['five_log']*TF2OverallAmpl*((pi*fInsJoin)^(5. /3.)) - TF2coeffs['four']*TF2OverallAmpl*((pi*fInsJoin)^(4. /3.)) - 2. *TF2coeffs['three']*TF2OverallAmpl*(pi*fInsJoin) - 3. *TF2coeffs['two']*TF2OverallAmpl*((pi*fInsJoin)^(2. /3.)) - 4. *TF2coeffs['one']*TF2OverallAmpl*((pi*fInsJoin)^(1. /3.)) - 5. *TF2coeffs['zero']*TF2OverallAmpl)*pi/(3. *((pi*fInsJoin)^(8. /3.)))
    DPhiIns =
        (
            2.0 * TF2coeffs.seven * TF2OverallAmpl * ((pi * fInsJoin)^(7.0 / 3.0)) +
            (
                TF2coeffs.six * TF2OverallAmpl +
                TF2coeffs.six_log * TF2OverallAmpl * (1.0 + log(pi * fInsJoin) / 3.0)
            ) * ((pi * fInsJoin)^(2.0)) +
            TF2coeffs.five_log * TF2OverallAmpl * ((pi * fInsJoin)^(5.0 / 3.0)) -
            TF2coeffs.four * TF2OverallAmpl * ((pi * fInsJoin)^(4.0 / 3.0)) -
            2.0 * TF2coeffs.three * TF2OverallAmpl * (pi * fInsJoin) -
            3.0 * TF2coeffs.two * TF2OverallAmpl * ((pi * fInsJoin)^(2.0 / 3.0)) -
            4.0 * TF2coeffs.one * TF2OverallAmpl * ((pi * fInsJoin)^(1.0 / 3.0)) -
            5.0 * TF2coeffs.zero * TF2OverallAmpl
        ) * pi / (3.0 * ((pi * fInsJoin)^(8.0 / 3.0))) +
        (
            sigma1 +
            sigma2 * (fInsJoin^(1.0 / 3.0)) +
            sigma3 * (fInsJoin^(2.0 / 3.0)) +
            sigma4 * fInsJoin
        ) * etaInv
    #DPhiIns = DPhiIns + (sigma1 + sigma2*(fInsJoin^(1. /3.)) + sigma3*(fInsJoin^(2. /3.)) + sigma4*fInsJoin)*etaInv
    # This is the first derivative of the Intermediate phase computed at fInsJoin
    DPhiInt = (beta1 + beta3 / (fInsJoin^4) + beta2 / fInsJoin) * etaInv
    C2Int = DPhiIns - DPhiInt

    # This is the inspiral phase computed at fInsJoin
    # PhiInsJoin = PhiInspcoeffs['initial_phasing'] + PhiInspcoeffs['two_thirds']*(fInsJoin^(2. /3.)) + PhiInspcoeffs['third']*(fInsJoin^(1. /3.)) + PhiInspcoeffs['third_log']*(fInsJoin^(1. /3.))*log(pi*fInsJoin)/3. + PhiInspcoeffs['log']*log(pi*fInsJoin)/3. + PhiInspcoeffs['min_third']*(fInsJoin^(-1. /3.)) + PhiInspcoeffs['min_two_thirds']*(fInsJoin^(-2. /3.)) + PhiInspcoeffs['min_one']/fInsJoin + PhiInspcoeffs['min_four_thirds']*(fInsJoin^(-4. /3.)) + PhiInspcoeffs['min_five_thirds']*(fInsJoin^(-5. /3.)) + (PhiInspcoeffs['one']*fInsJoin + PhiInspcoeffs['four_thirds']*(fInsJoin^(4. /3.)) + PhiInspcoeffs['five_thirds']*(fInsJoin^(5. /3.)) + PhiInspcoeffs['two']*fInsJoin*fInsJoin)/eta
    PhiInsJoin =
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (fInsJoin^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (fInsJoin^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (fInsJoin^(1.0 / 3.0)) * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.log * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.min_third * (fInsJoin^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (fInsJoin^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / fInsJoin +
        PhiInspcoeffs.min_four_thirds * (fInsJoin^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (fInsJoin^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * fInsJoin +
            PhiInspcoeffs.four_thirds * (fInsJoin^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (fInsJoin^(5.0 / 3.0)) +
            PhiInspcoeffs.two * fInsJoin * fInsJoin
        ) * etaInv
    # This is the Intermediate phase computed at fInsJoin
    PhiIntJoin =
        beta1 * fInsJoin - beta3 / (3.0 * fInsJoin * fInsJoin * fInsJoin) +
        beta2 * log(fInsJoin)

    C1Int = PhiInsJoin - PhiIntJoin * etaInv - C2Int * fInsJoin

    # Now the same for Intermediate - Merger-Ringdown: we also need a temporary Intermediate Phase function
    PhiIntTempVal =
        (beta1 * fMRDJoin - beta3 / (3.0 * fMRDJoin^3) + beta2 * log(fMRDJoin)) * etaInv +
        C1Int +
        C2Int * fMRDJoin
    DPhiIntTempVal = C2Int + (beta1 + beta3 / (fMRDJoin^4) + beta2 / fMRDJoin) * etaInv
    DPhiMRDVal =
        (
            alpha1 +
            alpha2 / (fMRDJoin^2) +
            alpha3 / (fMRDJoin^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fMRDJoin - alpha5 * fring) * (fMRDJoin - alpha5 * fring) / (fdamp^2)
                )
            )
        ) * etaInv
    PhiMRJoinTemp =
        -(alpha2 / fMRDJoin) +
        (4.0 / 3.0) * (alpha3 * (fMRDJoin^(0.75))) +
        alpha1 * fMRDJoin +
        alpha4 * atan((fMRDJoin - alpha5 * fring) / fdamp)

    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRD * fMRDJoin
    # Time shift so that peak amplitude is approximately at t=0
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi^2
        ) * xi
    gamma3 =
        1.3081615607036106 - 0.005537729694807678 * eta +
        (
            -0.06782917938621007 - 0.6689834970767117 * eta +
            3.403147966134083 * eta2 +
            (-0.05296577374411866 - 0.9923793203111362 * eta + 4.820681208409587 * eta2) *
            xi +
            (
                -0.006134139870393713 - 0.38429253308696365 * eta +
                1.7561754421985984 * eta2
            ) * xi^2
        ) * xi
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs.(fring - (fdamp * gamma3) / gamma2),
        abs.(fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2^2)) * gamma3) / gamma2),
    )
    t0 =
        (
            alpha1 +
            alpha2 / (fpeak^2) +
            alpha3 / (fpeak^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fpeak - alpha5 * fring) * (fpeak - alpha5 * fring) / (fdamp * fdamp)
                )
            )
        ) * etaInv

    # LAL sets fRef as the minimum frequency, do the same
    fRef = fgrid[1]   #minimum(fgrid)

    phiRef = ifelse(
        fRef < fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fRef^(2.0 / 3.0) +
        PhiInspcoeffs.third * fRef^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fRef^(1.0 / 3.0) * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.log * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.min_third * fRef^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fRef^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fRef +
        PhiInspcoeffs.min_four_thirds * fRef^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fRef^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fRef +
            PhiInspcoeffs.four_thirds * fRef^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fRef^(5.0 / 3.0) +
            PhiInspcoeffs.two * fRef * fRef
        ) * etaInv,
        ifelse(
            fRef < fMRDJoin,
            (beta1 * fRef - beta3 / (3.0 * fRef * fRef * fRef) + beta2 * log(fRef)) *
            etaInv +
            C1Int +
            C2Int * fRef,
            ifelse(
                fRef < fcutPar,
                (
                    -(alpha2 / fRef) +
                    (4.0 / 3.0) * (alpha3 * (fRef^(3.0 / 4.0))) +
                    alpha1 * fRef +
                    alpha4 * atan((fRef - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fRef,
                0.0,
            ),
        ),
    )
    
    phis = @. ifelse.(
        fgrid .< fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fgrid^(2.0 / 3.0) +
        PhiInspcoeffs.third * fgrid^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fgrid^(1.0 / 3.0) * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.log * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.min_third * fgrid^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fgrid^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fgrid +
        PhiInspcoeffs.min_four_thirds * fgrid^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fgrid^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fgrid +
            PhiInspcoeffs.four_thirds * fgrid^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fgrid^(5.0 / 3.0) +
            PhiInspcoeffs.two * fgrid * fgrid
        ) * etaInv,
        ifelse(
            fgrid < fMRDJoin,
            (beta1 * fgrid - beta3 / (3.0 * fgrid * fgrid * fgrid) + beta2 * log(fgrid)) * etaInv +
            C1Int +
            C2Int * fgrid,
            ifelse(
                fgrid < fcutPar,
                (
                    -(alpha2 / fgrid) +
                    (4.0 / 3.0) * (alpha3 * (fgrid^(3.0 / 4.0))) +
                    alpha1 * fgrid +
                    alpha4 * atan((fgrid - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fgrid,
                0.0,
            ),
        ),
    )

    # Add the tidal contribution to the phase, as in arXiv:1905.06011
    # Compute the tidal coupling constant, arXiv:1905.06011 eq. (8) using Lambda = 2/3 k_2/C^5 (eq. (10))

    kappa2T = (3.0/13.0) * ((1.0 + 12.0*m2ByM/m1ByM)*(m1ByM^5)*Lambda1 + (1.0 + 12.0*m1ByM/m2ByM)*(m2ByM^5)*Lambda2)
    
    c_Newt   = 2.4375
    n_1      = -12.615214237993088
    n_3over2 =  19.0537346970349
    n_2      = -21.166863146081035
    n_5over2 =  90.55082156324926
    n_3      = -60.25357801943598
    d_1      = -15.11120782773667
    d_3over2 =  22.195327350624694
    d_2      =   8.064109635305156

    numTidal = @.  1.0 + (n_1 * ((pi*fgrid)^(2. /3.))) + (n_3over2 * pi*fgrid) + (n_2 * ((pi*fgrid)^(4. /3.))) + (n_5over2 * ((pi*fgrid)^(5. /3.))) + (n_3 * (pi*fgrid)^2)
    denTidal = @. 1.0 + (d_1 * ((pi*fgrid)^(2. /3.))) + (d_3over2 * pi*fgrid) + (d_2 * ((pi*fgrid)^(4. /3.)))
    
    tidal_phase = @. - kappa2T * c_Newt / (m1ByM * m2ByM) * ((pi*fgrid)^(5. /3.)) * numTidal / denTidal
    
    # In the NRTidalv2 extension also 3.5PN spin-squared and 3.5PN spin-cubed terms are included, see eq. (27) of arXiv:1905.06011
    # This is needed to account for spin-induced quadrupole moments
    # Compute octupole moment and emove -1 to account for BBH baseline
    
    OctMon1 = - 1. + exp(0.003131 + 2.071 * log(QuadMon1)  - 0.7152 * log(QuadMon1)^2 + 0.2458 * log(QuadMon1)^3 - 0.03309 * log(QuadMon1)^4)
    OctMon2 = - 1. + exp(0.003131 + 2.071 * log(QuadMon2)  - 0.7152 * log(QuadMon2)^2 + 0.2458 * log(QuadMon2)^3 - 0.03309 * log(QuadMon2)^4)

    SS_3p5PN  = - 400. *pi*(QuadMon1-1.)*chi12*m1ByM^2 - 400. *pi*(QuadMon2-1.)*chi22*m2ByM^2
    SSS_3p5PN = 10. *((m1ByM^2+308. /3. *m1ByM)*chi1+(m2ByM^2-89. /3. *m2ByM)*chi2)*(QuadMon1-1.)*m1ByM^2*chi12 + 10. *((m2ByM^2+308. /3. *m2ByM)*chi2+(m1ByM^2-89. /3. *m1ByM)*chi1)*(QuadMon2-1.)*m2ByM^2*chi22 - 440. *OctMon1*m1ByM^3*chi12*chi1 - 440. *OctMon2*m2ByM^3*chi22*chi2

    return @. phis + ifelse(fgrid < fcutPar, - t0*(fgrid - fRef) - phiRef + tidal_phase + (SS_3p5PN + SSS_3p5PN)*TF2OverallAmpl*((pi*fgrid)^(2. /3.)), 0.)
end  

"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomD_NRTidal(), f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `dL`: Luminosity distance to the source, in Gpc.
-  `Lambda1`: Dimensionless tidal deformability of the first NS.
-  `Lambda2`: Dimensionless tidal deformability of the second NS.
#### Optional arguments:
-  `fInsJoin`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The amplitude is dimensionless.

#### Example:
```julia
    mc = 3.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    Lambda1 = 1000.
    Lambda2 = 2000.
    fcut = _fcut(PhenomD_NRTidal(), mc, eta, Lambda1, Lambda2)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomD_NRTidal(), f, mc, eta, chi1, chi2, dL, Lambda1, Lambda2)
```
"""
function Ampl(model::PhenomD_NRTidal,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    # # print that if Lambda1 and or Lambda2 are zero one should resolve to other waveforms
    # if Lambda1 == 0. && Lambda2 ==0.
    #     println("Lambda1 and Lambda2 are zero, please use other waveforms, valid for BBH")
    # elseif Lambda1 == 0. || Lambda2 == 0.
    #     println("Lambda1 or Lambda2 is zero, please use PhenomNSBH")
    # end

    ampl_notidal = Ampl(PhenomD(),
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = fcutPar,
    fInsJoin_Ampl = fInsJoin_Ampl
    )  
    
    M = mc/eta^(3. /5.)
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.) 
    fgrid = M * GMsun_over_c3 .* f
    Overallamp =
    2.0 * sqrt(5.0 / (64.0 * pi)) * M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL

    # Now add the tidal amplitude as in arXiv:1905.06011
    # Compute the tidal coupling constant, arXiv:1905.06011 eq. (8) using Lambda = 2/3 k_2/C^5 (eq. (10))
    Xa = 0.5 * (1.0 + Seta)
    Xb = 0.5 * (1.0 - Seta)
    kappa2T = (3.0/13.0) * ((1.0 + 12.0*Xb/Xa)*(Xa^5)*Lambda1 + (1.0 + 12.0*Xa/Xb)*(Xb^5)*Lambda2)
    
    # Now compute the amplitude modification as in arXiv:1905.06011 eq. (24)
    xTidal = @. (pi * fgrid)^(2. /3.)
    n1T    = 4.157407407407407
    n289T  = 2519.111111111111
    dTidal = 13477.8073677
    polyTidal = @. (1.0 + n1T*xTidal + n289T*(xTidal^(2.89)))/(1. +dTidal*(xTidal^4))
    ampTidal = @. -9.0*kappa2T*(xTidal^3.25)*polyTidal
    
    # # Compute the dimensionless merger frequency (Mf) for the Planck taper filtering
    # a_0 = 0.3586
    # n_1 = 3.35411203e-2
    # n_2 = 4.31460284e-5
    # d_1 = 7.54224145e-2
    # d_2 = 2.23626859e-4

    
    # numPT = 1.0 + n_1*kappa2T + n_2*kappa2T^2
    # denPT = 1.0 + d_1*kappa2T + d_2*kappa2T^2
    # Q_0 = a_0 / sqrt(q)
    # f_merger = Q_0 * (numPT / denPT) / (2. *pi)
    f_end_taper = _fcut(model, mc, eta, Lambda1, Lambda2) * GMsun_over_c3 * M
    f_merger = f_end_taper / 1.2
    # The derivative of the Planck taper filter can return NaN in some points because of numerical issues, we declare it explicitly to avoid the issue # ANDREA check
    planck_taper = @. ifelse(fgrid < f_merger, 1., ifelse(fgrid > f_end_taper, 0., 1. - 1. /(exp((f_end_taper - f_merger)/(fgrid - f_merger) + (f_end_taper - f_merger)/(fgrid - f_end_taper)) + 1.)))

    return  @. (ampl_notidal + 2*sqrt(pi/5.)*ampTidal*Overallamp)*planck_taper
end 



"""
Compute the cut frequency of the waveform as a function of the events parameters, in `Hz`.
Valid for PhenomD_NRTidal.
We cut the waveform slightly before the end of the Planck taper filter, for numerical stability.

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `Lambda1`: Dimensionless tidal deformability of the first NS.
-  `Lambda2`: Dimensionless tidal deformability of the second NS.

#### Return:
-  (float) Cut frequency of the waveform, in Hz.

"""
  
function _fcut(model::PhenomD_NRTidal, mc, eta, Lambda1, Lambda2; GMsun_over_c3=uc.GMsun_over_c3)

    M = mc/eta^(3. /5.)
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    Xa = 0.5 * (1.0 + Seta)
    Xb = 0.5 * (1.0 - Seta)
    kappa2T = (3.0/13.0) * ((1.0 + 12.0*Xb/Xa)*(Xa^5)*Lambda1 + (1.0 + 12.0*Xa/Xb)*(Xb^5)*Lambda2)
    
    # Compute the dimensionless merger frequency (Mf) for the Planck taper filtering
    a_0 = 0.3586
    n_1 = 3.35411203e-2
    n_2 = 4.31460284e-5
    d_1 = 7.54224145e-2
    d_2 = 2.23626859e-4
    
    numPT = 1.0 + n_1*kappa2T + n_2*kappa2T^2
    denPT = 1.0 + d_1*kappa2T + d_2*kappa2T^2
    Q_0 = a_0 / sqrt(q)
    f_merger = Q_0 * (numPT / denPT) / (2. *pi)
    # Terminate the waveform at 1.2 times the merger frequency
    f_end_taper = 1.2*f_merger

    return f_end_taper/(M*GMsun_over_c3)
end



############################################################################################################
#   IMRPHENOM  HM 
############################################################################################################
"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomHM(), f, mc, eta, chi1, chi2)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts an array of frequencies.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Optional arguments:
-  `fInsJoin`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns a structure containg the 6 modes of the phase. The phase is given in radians.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomHM(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Phi(PhenomHM(), f, mc, eta, chi1, chi2)
```
"""

function Phi(model::PhenomHM,
    f,
    mc,
    eta,
    chi1,
    chi2;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)


    complShiftm = [0.0, pi * 0.5, 0.0, -pi * 0.5, pi, pi * 0.5, 0.0]

    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta
    pi2 = pi * pi

    QuadMon1, QuadMon2 = 1.0, 1.0

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    q = 0.5 * (1.0 + Seta - 2.0 * eta) * etaInv
    chi_s2, chi_a2 = chi_s * chi_s, chi_a * chi_a
    chi1dotchi2 = chi1 * chi2
    chi_sdotchi_a = chi_s * chi_a
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # This is MfRef, needed to recover LAL, which sets fRef to f_min if fRef=0
    fRef = fgrid[1] #minimum(fgrid)
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin, radiated energy and mass
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    finMass = 1.0 - Erad

    fring, fdamp = _RDfreqCalc(model, finMass, aeff, 2, 2)

    # Compute sigma coefficients appearing in arXiv:1508.07253 eq. (28)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    sigma1 =
        2096.551999295543 +
        1463.7493168261553 * eta +
        (
            1312.5493286098522 + 18307.330017082117 * eta - 43534.1440746107 * eta2 +
            (-833.2889543511114 + 32047.31997183187 * eta - 108609.45037520859 * eta2) *
            xi +
            (452.25136398112204 + 8353.439546391714 * eta - 44531.3250037322 * eta2) * xi2
        ) * xi
    sigma2 =
        -10114.056472621156 - 44631.01109458185 * eta +
        (
            -6541.308761668722 - 266959.23419307504 * eta +
            686328.3229317984 * eta2 +
            (3405.6372187679685 - 437507.7208209015 * eta + 1.6318171307344697e6 * eta2) *
            xi +
            (-7462.648563007646 - 114585.25177153319 * eta + 674402.4689098676 * eta2) * xi2
        ) * xi
    sigma3 =
        22933.658273436497 +
        230960.00814979506 * eta +
        (
            14961.083974183695 + 1.1940181342318142e6 * eta - 3.1042239693052764e6 * eta2 +
            (-3038.166617199259 + 1.8720322849093592e6 * eta - 7.309145012085539e6 * eta2) *
            xi +
            (42738.22871475411 + 467502.018616601 * eta - 3.064853498512499e6 * eta2) * xi2
        ) * xi
    sigma4 =
        -14621.71522218357 - 377812.8579387104 * eta +
        (
            -9608.682631509726 - 1.7108925257214056e6 * eta +
            4.332924601416521e6 * eta2 +
            (
                -22366.683262266528 - 2.5019716386377467e6 * eta +
                1.0274495902259542e7 * eta2
            ) * xi +
            (-85360.30079034246 - 570025.3441737515 * eta + 4.396844346849777e6 * eta2) *
            xi2
        ) * xi

    # Compute beta coefficients appearing in arXiv:1508.07253 eq. (16)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    beta1 =
        97.89747327985583 - 42.659730877489224 * eta +
        (
            153.48421037904913 - 1417.0620760768954 * eta +
            2752.8614143665027 * eta2 +
            (138.7406469558649 - 1433.6585075135881 * eta + 2857.7418952430758 * eta2) *
            xi +
            (41.025109467376126 - 423.680737974639 * eta + 850.3594335657173 * eta2) * xi2
        ) * xi
    beta2 =
        -3.282701958759534 - 9.051384468245866 * eta +
        (
            -12.415449742258042 + 55.4716447709787 * eta - 106.05109938966335 * eta2 +
            (-11.953044553690658 + 76.80704618365418 * eta - 155.33172948098394 * eta2) *
            xi +
            (-3.4129261592393263 + 25.572377569952536 * eta - 54.408036707740465 * eta2) *
            xi2
        ) * xi
    beta3 =
        -0.000025156429818799565 +
        0.000019750256942201327 * eta +
        (
            -0.000018370671469295915 +
            0.000021886317041311973 * eta +
            0.00008250240316860033 * eta2 +
            (
                7.157371250566708e-6 - 0.000055780000112270685 * eta +
                0.00019142082884072178 * eta2
            ) * xi +
            (
                5.447166261464217e-6 - 0.00003220610095021982 * eta +
                0.00007974016714984341 * eta2
            ) * xi2
        ) * xi

    # Compute alpha coefficients appearing in arXiv:1508.07253 eq. (14)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    alpha1 =
        43.31514709695348 +
        638.6332679188081 * eta +
        (
            -32.85768747216059 + 2415.8938269370315 * eta - 5766.875169379177 * eta2 +
            (-61.85459307173841 + 2953.967762459948 * eta - 8986.29057591497 * eta2) * xi +
            (-21.571435779762044 + 981.2158224673428 * eta - 3239.5664895930286 * eta2) *
            xi2
        ) * xi
    alpha2 =
        -0.07020209449091723 - 0.16269798450687084 * eta +
        (
            -0.1872514685185499 + 1.138313650449945 * eta - 2.8334196304430046 * eta2 +
            (-0.17137955686840617 + 1.7197549338119527 * eta - 4.539717148261272 * eta2) *
            xi +
            (-0.049983437357548705 + 0.6062072055948309 * eta - 1.682769616644546 * eta2) *
            xi2
        ) * xi
    alpha3 =
        9.5988072383479 - 397.05438595557433 * eta +
        (
            16.202126189517813 - 1574.8286986717037 * eta +
            3600.3410843831093 * eta2 +
            (27.092429659075467 - 1786.482357315139 * eta + 5152.919378666511 * eta2) * xi +
            (11.175710130033895 - 577.7999423177481 * eta + 1808.730762932043 * eta2) * xi2
        ) * xi
    alpha4 =
        -0.02989487384493607 +
        1.4022106448583738 * eta +
        (
            -0.07356049468633846 +
            0.8337006542278661 * eta +
            0.2240008282397391 * eta2 +
            (-0.055202870001177226 + 0.5667186343606578 * eta + 0.7186931973380503 * eta2) *
            xi +
            (
                -0.015507437354325743 +
                0.15750322779277187 * eta +
                0.21076815715176228 * eta2
            ) * xi2
        ) * xi
    alpha5 =
        0.9974408278363099 - 0.007884449714907203 * eta +
        (
            -0.059046901195591035 + 1.3958712396764088 * eta - 4.516631601676276 * eta2 +
            (-0.05585343136869692 + 1.7516580039343603 * eta - 5.990208965347804 * eta2) *
            xi +
            (-0.017945336522161195 + 0.5965097794825992 * eta - 2.0608879367971804 * eta2) *
            xi2
        ) * xi

    # Compute the TF2 phase coefficients and put them in a dictionary (spin effects are included up to 3.5PN)

    TF2OverallAmpl = 3 / (128.0 * eta)
    TF2_5coeff_tmp =
        38645.0 * pi / 756.0 - 65.0 * pi * eta / 9.0 - (
            (732985.0 / 2268.0 - 24260.0 * eta / 81.0 - 340.0 * eta2 / 9.0) * chi_s +
            (732985.0 / 2268.0 + 140.0 * eta / 9.0) * Seta * chi_a
        ) #variable to be used later
    TF2_6coeff_tmp =
        11583.231236531 / 4.694215680 - 640.0 / 3.0 * pi2 -
        684.8 / 2.1 * MathConstants.eulergamma +
        eta * (-15737.765635 / 3.048192 + 225.5 / 1.2 * pi2) +
        eta2 * 76.055 / 1.728 - eta2 * eta * 127.825 / 1.296 - log(4.0) * 684.8 / 2.1 +
        pi * chi1 * m1ByM * (1490.0 / 3.0 + m1ByM * 260.0) +
        pi * chi2 * m2ByM * (1490.0 / 3.0 + m2ByM * 260.0) +
        (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) *
        m1ByM^2 *
        QuadMon1 *
        chi12 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 ) *
        m1ByM^2 *
        chi12 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) *
        m2ByM^2 *
        QuadMon2 *
        chi22 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 ) *
        m2ByM^2 *
        chi22

    TF2coeffs = TF2coeffsStructure(
        1.0,
        0.0,
        3715.0 / 756.0 + (55.0 * eta) / 9.0,
        -16.0 * pi +
        (113.0 * Seta * chi_a) / 3.0 +
        (113.0 / 3.0 - (76.0 * eta) / 3.0) * chi_s,
        5.0 * (3058.673 / 7.056 + 5429.0 / 7.0 * eta + 617.0 * eta2) / 72.0 +
        247.0 / 4.8 * eta * chi1dotchi2 - 721.0 / 4.8 * eta * chi1dotchi2 +
        (-720.0 / 9.6 * QuadMon1 + 1.0 / 9.6) * m1ByM^2  * chi12 +
        (-720.0 / 9.6 * QuadMon2 + 1.0 / 9.6) * m2ByM^2  * chi22 +
        (240.0 / 9.6 * QuadMon1 - 7.0 / 9.6) * m1ByM^2  * chi12 +
        (240.0 / 9.6 * QuadMon2 - 7.0 / 9.6) * m2ByM^2  * chi22,
        TF2_5coeff_tmp,
        TF2_5coeff_tmp * 3.0,
        TF2_6coeff_tmp - (
            (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 )
            ) *
            m1ByM *
            m1ByM *
            chi12 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 )
            ) *
            m2ByM *
            m2ByM *
            chi22
        ),
        -6848.0 / 21.0,
        77096675.0 * pi / 254016.0 + 378515.0 * pi * eta / 1512.0 -
        74045.0 * pi * eta2 / 756.0 +
        (
            -25150083775.0 / 3048192.0 + 10566655595.0 * eta / 762048.0 -
            1042165.0 * eta2 / 3024.0 + 5345.0 * eta2 * eta / 36.0
        ) * chi_s +
        Seta * (
            (
                -25150083775.0 / 3048192.0 + 26804935.0 * eta / 6048.0 -
                1985.0 * eta2 / 48.0
            ) * chi_a
        ),
    )


    PhiInspcoeffs = PhiInspcoeffsStructure(
        TF2coeffs.five * TF2OverallAmpl,
        TF2coeffs.seven * TF2OverallAmpl * (pi^(2.0 / 3.0)),
        TF2coeffs.six * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.six_log * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.five_log * TF2OverallAmpl,
        TF2coeffs.four * TF2OverallAmpl * (pi^(-1.0 / 3.0)),
        TF2coeffs.three * TF2OverallAmpl * (pi^(-2.0 / 3.0)),
        TF2coeffs.two * TF2OverallAmpl / pi,
        TF2coeffs.one * TF2OverallAmpl * (pi^(-4.0 / 3.0)),
        TF2coeffs.zero * TF2OverallAmpl * (pi^(-5.0 / 3.0)),
        sigma1,
        sigma2 * 0.75,
        sigma3 * 0.6,
        sigma4 * 0.5,
    )

    #Now compute the coefficients to align the three parts

    fInsJoin = fInsJoin_PHI
    fMRDJoin = 0.5 * fring

    # First the Inspiral - Intermediate: we compute C1Int and C2Int coeffs
    # Equations to solve for to get C(1) continuous join
    # PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    # Joining at fInsJoin
    # PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    # PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    # This is the first derivative wrt f of the inspiral phase computed at fInsJoin, first add the PN contribution and then the higher order calibrated terms
    DPhiIns =
        (
            2.0 * TF2coeffs.seven * TF2OverallAmpl * ((pi * fInsJoin)^(7.0 / 3.0)) +
            (
                TF2coeffs.six * TF2OverallAmpl +
                TF2coeffs.six_log * TF2OverallAmpl * (1.0 + log(pi * fInsJoin) / 3.0)
            ) * ((pi * fInsJoin)^(2.0)) +
            TF2coeffs.five_log * TF2OverallAmpl * ((pi * fInsJoin)^(5.0 / 3.0)) -
            TF2coeffs.four * TF2OverallAmpl * ((pi * fInsJoin)^(4.0 / 3.0)) -
            2.0 * TF2coeffs.three * TF2OverallAmpl * (pi * fInsJoin) -
            3.0 * TF2coeffs.two * TF2OverallAmpl * ((pi * fInsJoin)^(2.0 / 3.0)) -
            4.0 * TF2coeffs.one * TF2OverallAmpl * ((pi * fInsJoin)^(1.0 / 3.0)) -
            5.0 * TF2coeffs.zero * TF2OverallAmpl
        ) * pi / (3.0 * ((pi * fInsJoin)^(8.0 / 3.0))) +
        (
            sigma1 +
            sigma2 * (fInsJoin^(1.0 / 3.0)) +
            sigma3 * (fInsJoin^(2.0 / 3.0)) +
            sigma4 * fInsJoin
        ) * etaInv
    # This is the first derivative of the Intermediate phase computed at fInsJoin
    DPhiInt = (beta1 + beta3 / (fInsJoin^4) + beta2 / fInsJoin) * etaInv

    C2Int = DPhiIns - DPhiInt

    # This is the inspiral phase computed at fInsJoin 
    PhiInsJoin =
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (fInsJoin^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (fInsJoin^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (fInsJoin^(1.0 / 3.0)) * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.log * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.min_third * (fInsJoin^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (fInsJoin^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / fInsJoin +
        PhiInspcoeffs.min_four_thirds * (fInsJoin^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (fInsJoin^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * fInsJoin +
            PhiInspcoeffs.four_thirds * (fInsJoin^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (fInsJoin^(5.0 / 3.0)) +
            PhiInspcoeffs.two * fInsJoin * fInsJoin
        ) * etaInv
    # This is the Intermediate phase computed at fInsJoin
    PhiIntJoin =
        beta1 * fInsJoin - beta3 / (3.0 * fInsJoin * fInsJoin * fInsJoin) +
        beta2 * log(fInsJoin)

    C1Int = PhiInsJoin - PhiIntJoin * etaInv - C2Int * fInsJoin

    # Now the same for Intermediate - Merger-Ringdown: we also need a temporary Intermediate Phase function
    PhiIntTempVal =
        (beta1 * fMRDJoin - beta3 / (3.0 * fMRDJoin^3) + beta2 * log(fMRDJoin)) * etaInv +
        C1Int +
        C2Int * fMRDJoin
    DPhiIntTempVal = C2Int + (beta1 + beta3 / (fMRDJoin^4) + beta2 / fMRDJoin) * etaInv
    DPhiMRDVal =
        (
            alpha1 +
            alpha2 / (fMRDJoin^2) +
            alpha3 / (fMRDJoin^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fMRDJoin - alpha5 * fring) * (fMRDJoin - alpha5 * fring) / (fdamp^2)
                )
            )
        ) * etaInv
    PhiMRJoinTemp =
        -(alpha2 / fMRDJoin) +
        (4.0 / 3.0) * (alpha3 * (fMRDJoin^(0.75))) +
        alpha1 * fMRDJoin +
        alpha4 * atan((fMRDJoin - alpha5 * fring) / fdamp)

    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRD * fMRDJoin

    # Time shift so that peak amplitude is approximately at t=0
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi^2
        ) * xi
    gamma3 =
        1.3081615607036106 - 0.005537729694807678 * eta +
        (
            -0.06782917938621007 - 0.6689834970767117 * eta +
            3.403147966134083 * eta2 +
            (-0.05296577374411866 - 0.9923793203111362 * eta + 4.820681208409587 * eta2) *
            xi +
            (
                -0.006134139870393713 - 0.38429253308696365 * eta +
                1.7561754421985984 * eta2
            ) * xi^2
        ) * xi
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs.(fring - (fdamp * gamma3) / gamma2),
        abs.(fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2^2)) * gamma3) / gamma2),
    )
    t0 =
        (
            alpha1 +
            alpha2 / (fpeak^2) +
            alpha3 / (fpeak^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fpeak - alpha5 * fring) * (fpeak - alpha5 * fring) / (fdamp * fdamp)
                )
            )
        ) * etaInv


    phiRef = _completePhase(
        fRef,
        C1MRD,
        C2MRD,
        1.0,
        1.0,
        C1Int,
        C2Int,
        PhiInspcoeffs,
        fInsJoin,
        alpha1,
        alpha2,
        alpha3,
        alpha4,
        alpha5,
        beta1,
        beta2,
        beta3,
        etaInv,
        fMRDJoin,
        fcutPar,
        fring,
        fdamp,
    )
    phi0 = 0.5 * phiRef

    # Now compute the other modes, they are 5, we loop
    tmpResults = []
    for ell in (2:4)
        for mm in (ell - 1, ell)
            # Compute ringdown and damping frequencies from fits for the various modes
            fringlm, fdamplm = _RDfreqCalc(model, finMass, aeff, ell, mm)
            Rholm, Taulm = fring / fringlm, fdamplm / fdamp

            # Rholm and Taulm only figure in the MRD part, the rest of the coefficients is the same, recompute only this
            DPhiMRDVal =
                (
                    alpha1 +
                    alpha2 / (fMRDJoin * fMRDJoin) +
                    alpha3 / (fMRDJoin^(1.0 / 4.0)) +
                    alpha4 / (
                        fdamp *
                        Taulm *
                        (
                            1.0 +
                            (fMRDJoin - alpha5 * fring) * (fMRDJoin - alpha5 * fring) /
                            (fdamp * Taulm * Rholm * fdamp * Taulm * Rholm)
                        )
                    )
                ) * etaInv
            PhiMRJoinTemp =
                -(alpha2 / fMRDJoin) +
                (4.0 / 3.0) * (alpha3 * (fMRDJoin^(3.0 / 4.0))) +
                alpha1 * fMRDJoin +
                alpha4 * Rholm * atan((fMRDJoin - alpha5 * fring) / (fdamp * Rholm * Taulm))
            C2MRDHM = DPhiIntTempVal - DPhiMRDVal
            C1MRDHM = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRDHM * fMRDJoin

            # Compute mapping coefficinets
            Map_fl = fInsJoin
            Map_fi = Map_fl / Rholm
            Map_fr = fringlm

            Map_ai, Map_bi = 2.0 / mm, 0.0

            Map_Trd = Map_fr * Rholm
            Map_Ti = 2.0 * Map_fi / mm
            Map_am = (Map_Trd - Map_Ti) / (Map_fr - Map_fi)
            Map_bm = Map_Ti - Map_fi * Map_am

            Map_ar, Map_br = Rholm, 0.0
            # Now compute the needed constants
            tmpMf = Map_am * Map_fi + Map_bm
            PhDBconst =
                1.0 ./ Map_am .* _completePhase(
                    tmpMf,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoin,
                    fcutPar,
                    fring,
                    fdamp,
                )

            tmpMf = Map_ar * Map_fr + Map_br
            PhDCconst =
                1.0 ./ Map_ar .* _completePhase(
                    tmpMf,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoin,
                    fcutPar,
                    fring,
                    fdamp,
                )

            tmpMf = Map_ai * Map_fi + Map_bi
            PhDBAterm =
                1.0 ./ Map_ai .* _completePhase(
                    tmpMf,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoin,
                    fcutPar,
                    fring,
                    fdamp,
                )

            tmpMf = Map_am * Map_fr + Map_bm
            tmpphaseC =
                -PhDBconst .+ PhDBAterm .+
                1.0 ./ Map_am .* _completePhase(
                    tmpMf,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoin,
                    fcutPar,
                    fring,
                    fdamp,
                )

            phitmp =
                ifelse.(
                    fgrid .< Map_fi,
                    1.0 ./ Map_ai .* _completePhase(
                        fgrid .* Map_ai .+ Map_bi,
                        C1MRDHM,
                        C2MRDHM,
                        Rholm,
                        Taulm,
                        C1Int,
                        C2Int,
                        PhiInspcoeffs,
                        fInsJoin,
                        alpha1,
                        alpha2,
                        alpha3,
                        alpha4,
                        alpha5,
                        beta1,
                        beta2,
                        beta3,
                        etaInv,
                        fMRDJoin,
                        fcutPar,
                        fring,
                        fdamp,
                    ),
                    ifelse.(
                        fgrid .< Map_fr,
                        -PhDBconst .+ PhDBAterm .+
                        1.0 ./ Map_am .* _completePhase(
                            fgrid .* Map_am .+ Map_bm,
                            C1MRDHM,
                            C2MRDHM,
                            Rholm,
                            Taulm,
                            C1Int,
                            C2Int,
                            PhiInspcoeffs,
                            fInsJoin,
                            alpha1,
                            alpha2,
                            alpha3,
                            alpha4,
                            alpha5,
                            beta1,
                            beta2,
                            beta3,
                            etaInv,
                            fMRDJoin,
                            fcutPar,
                            fring,
                            fdamp,
                        ),
                        -PhDCconst .+ tmpphaseC .+
                        1.0 ./ Map_ar .* _completePhase(
                            fgrid .* Map_ar .+ Map_br,
                            C1MRDHM,
                            C2MRDHM,
                            Rholm,
                            Taulm,
                            C1Int,
                            C2Int,
                            PhiInspcoeffs,
                            fInsJoin,
                            alpha1,
                            alpha2,
                            alpha3,
                            alpha4,
                            alpha5,
                            beta1,
                            beta2,
                            beta3,
                            etaInv,
                            fMRDJoin,
                            fcutPar,
                            fring,
                            fdamp,
                        ),
                    ),
                )
            push!(
                tmpResults,
                @. phitmp - t0 * (fgrid - fRef) - mm * phi0 + complShiftm[mm+1]
            )
        end
    end
    phislmp = PhislmpStructure(
        tmpResults[1],
        tmpResults[2],
        tmpResults[3],
        tmpResults[4],
        tmpResults[5],
        tmpResults[6],
    )  # probably improvable this step
    # the structure contains 6 terms each and each of them is an array of the same length as fgrid
    return phislmp
end


"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomHM(), f, mc, eta, chi1, chi2, dL)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts an array of frequencies.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `dL`: Luminosity distance to the binary, in Gpc.
#### Optional arguments:
-  `fInsJoin_Ampl`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns a structure containg the 6 modes of the amplitude. 
           To access the amplitude of the (2,2) mode, for example, you can use `Ampl(PhenomHM(), f, mc, eta, dL, chi1, chi2).two_two`. See also the julia function `fieldnames` to see all the available fields.
The amplitude is dimensionless. 
#### Example:
```julia
    mc = 30.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomHM(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomD(), f, mc, eta, chi1, chi2, dL)
```
"""
function Ampl(model::PhenomHM,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)


    # Useful quantities
    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta
    pi2 = pi * pi

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin and radiated energy
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    finMass = 1.0 - Erad

    fring, fdamp = _RDfreqCalc(model, finMass, aeff, 2, 2)

    # Compute coefficients gamma appearing in arXiv:1508.07253 eq. (19), the numerical coefficients are in Tab. 5
    gamma1 =
        0.006927402739328343 +
        0.03020474290328911 * eta +
        (
            0.006308024337706171 - 0.12074130661131138 * eta +
            0.26271598905781324 * eta2 +
            (
                0.0034151773647198794 - 0.10779338611188374 * eta +
                0.27098966966891747 * eta2
            ) * xi +
            (
                0.0007374185938559283 - 0.02749621038376281 * eta +
                0.0733150789135702 * eta2
            ) * xi2
        ) * xi
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi2
        ) * xi
    gamma3 =
        1.3081615607036106 - 0.005537729694807678 * eta +
        (
            -0.06782917938621007 - 0.6689834970767117 * eta +
            3.403147966134083 * eta2 +
            (-0.05296577374411866 - 0.9923793203111362 * eta + 4.820681208409587 * eta2) *
            xi +
            (
                -0.006134139870393713 - 0.38429253308696365 * eta +
                1.7561754421985984 * eta2
            ) * xi2
        ) * xi
    # Compute fpeak, from arXiv:1508.07253 eq. (20), we remove the square root term in case it is complex
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs(fring - (fdamp * gamma3) / gamma2),
        fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2 * gamma2)) * gamma3) / gamma2,
    )
    # Compute coefficients rho appearing in arXiv:1508.07253 eq. (30), the numerical coefficients are in Tab. 5
    rho1 =
        3931.8979897196696 - 17395.758706812805 * eta +
        (
            3132.375545898835 + 343965.86092361377 * eta - 1.2162565819981997e6 * eta2 +
            (-70698.00600428853 + 1.383907177859705e6 * eta - 3.9662761890979446e6 * eta2) *
            xi +
            (-60017.52423652596 + 803515.1181825735 * eta - 2.091710365941658e6 * eta2) *
            xi2
        ) * xi
    rho2 =
        -40105.47653771657 +
        112253.0169706701 * eta +
        (
            23561.696065836168 - 3.476180699403351e6 * eta +
            1.137593670849482e7 * eta2 +
            (754313.1127166454 - 1.308476044625268e7 * eta + 3.6444584853928134e7 * eta2) *
            xi +
            (596226.612472288 - 7.4277901143564405e6 * eta + 1.8928977514040343e7 * eta2) *
            xi2
        ) * xi
    rho3 =
        83208.35471266537 - 191237.7264145924 * eta +
        (
            -210916.2454782992 + 8.71797508352568e6 * eta - 2.6914942420669552e7 * eta2 +
            (
                -1.9889806527362722e6 + 3.0888029960154563e7 * eta -
                8.390870279256162e7 * eta2
            ) * xi +
            (
                -1.4535031953446497e6 + 1.7063528990822166e7 * eta -
                4.2748659731120914e7 * eta2
            ) * xi2
        ) * xi
    # Compute coefficients delta appearing in arXiv:1508.07253 eq. (21)
    f1Interm = fInsJoin_Ampl
    f3Interm = fpeak
    dfInterm = 0.5 * (f3Interm - f1Interm)
    f2Interm = f1Interm + dfInterm
    # First write the inspiral coefficients, we put them in a dictionary and label with the power in front of which they appear
    amp0 = sqrt(2.0 * eta / 3.0) * (pi^(-1.0 / 6.0))
    Acoeffs = AcoeffsStructure(
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / 672.0,
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0,
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) * (pi^(4.0 / 3.0))
        ) / 8.128512e6,
        (
            (pi^(5.0 / 3.0)) * (
                chi2 * (
                    -285197.0 * (-1.0 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1.0 + 4.0 * eta) * pi
            )
        ) / 32256.0,
        -(
            (pi^2.0) * (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368.0 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320.0 * pi2
                )
            )
        ) / 6.0085960704e10,
        rho1,
        rho2,
        rho3,
    )

    # v1 is the inspiral model evaluated at f1Interm
    v1 =
        1.0 +
        (f1Interm^(2.0 / 3.0)) * Acoeffs.two_thirds +
        (f1Interm^(4.0 / 3.0)) * Acoeffs.four_thirds +
        (f1Interm^(5.0 / 3.0)) * Acoeffs.five_thirds +
        (f1Interm^(7.0 / 3.0)) * Acoeffs.seven_thirds +
        (f1Interm^(8.0 / 3.0)) * Acoeffs.eight_thirds +
        f1Interm *
        (Acoeffs.one + f1Interm * Acoeffs.two + f1Interm * f1Interm * Acoeffs.three)
    # d1 is the derivative of the inspiral model evaluated at f1
    d1 =
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / (1008.0 * (f1Interm^(1.0 / 3.0))) +
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0 +
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) *
            (f1Interm^(1.0 / 3.0)) *
            (pi^(4.0 / 3.0))
        ) / 6.096384e6 +
        (
            5.0 *
            (f1Interm^(2.0 / 3.0)) *
            (pi^(5.0 / 3.0)) *
            (
                chi2 * (
                    -285197.0 * (-1 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1 + 4 * eta) * pi
            )
        ) / 96768.0 -
        (
            f1Interm *
            pi2 *
            (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320 * pi2
                )
            )
        ) / 3.0042980352e10 +
        (7.0 / 3.0) * (f1Interm^(4.0 / 3.0)) * rho1 +
        (8.0 / 3.0) * (f1Interm^(5.0 / 3.0)) * rho2 +
        3.0 * (f1Interm * f1Interm) * rho3
    # v3 is the merger-ringdown model (eq. (19) of arXiv:1508.07253) evaluated at f3
    v3 =
        exp(-(f3Interm - fring) * gamma2 / (fdamp * gamma3)) * (fdamp * gamma3 * gamma1) /
        ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
    # d2 is the derivative of the merger-ringdown model evaluated at f3
    d2 =
        (
            (-2.0 * fdamp * (f3Interm - fring) * gamma3 * gamma1) /
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3) -
            (gamma2 * gamma1)
        ) / (
            exp((f3Interm - fring) * gamma2 / (fdamp * gamma3)) *
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
        )
    # v2 is the value of the amplitude evaluated at f2. They come from the fit of the collocation points in the intermediate region
    v2 =
        0.8149838730507785 +
        2.5747553517454658 * eta +
        (
            1.1610198035496786 - 2.3627771785551537 * eta +
            6.771038707057573 * eta2 +
            (0.7570782938606834 - 2.7256896890432474 * eta + 7.1140380397149965 * eta2) *
            xi +
            (0.1766934149293479 - 0.7978690983168183 * eta + 2.1162391502005153 * eta2) *
            xi2
        ) * xi
    # Now some definitions to speed up
    f1 = f1Interm
    f2 = f2Interm
    f3 = f3Interm
    f12 = f1Interm * f1Interm
    f13 = f1Interm * f12
    f14 = f1Interm * f13
    f15 = f1Interm * f14
    f22 = f2Interm * f2Interm
    f23 = f2Interm * f22
    f24 = f2Interm * f23
    f32 = f3Interm * f3Interm
    f33 = f3Interm * f32
    f34 = f3Interm * f33
    f35 = f3Interm * f34
    # Finally conpute the deltas
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2.0 * d1 * f13 * f22 * f33 - 2.0 * d2 * f13 * f22 * f33 -
            d1 * f12 * f23 * f33 + d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 -
            d1 * f13 * f2 * f34 - d1 * f12 * f22 * f34 +
            2.0 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4.0 * f12 * f23 * f32 * v1 -
            3.0 * f1 * f24 * f32 * v1 - 8.0 * f12 * f22 * f33 * v1 +
            4.0 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4.0 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2.0 * f23 * f34 * v1 - 2.0 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3.0 * f14 * f33 * v2 -
            3.0 * f13 * f34 * v2 + f12 * f35 * v2 - f15 * f22 * v3 +
            2.0 * f14 * f23 * v3 - f13 * f24 * v3 + 2.0 * f15 * f2 * f3 * v3 -
            f14 * f22 * f3 * v3 - 4.0 * f13 * f23 * f3 * v3 +
            3.0 * f12 * f24 * f3 * v3 - 4.0 * f14 * f2 * f32 * v3 +
            8.0 * f13 * f22 * f32 * v3 - 4.0 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2 * d1 * f13 * f22 * f33 - 2 * d2 * f13 * f22 * f33 - d1 * f12 * f23 * f33 +
            d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 - d1 * f13 * f2 * f34 -
            d1 * f12 * f22 * f34 +
            2 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4 * f12 * f23 * f32 * v1 -
            3 * f1 * f24 * f32 * v1 - 8 * f12 * f22 * f33 * v1 +
            4 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2 * f23 * f34 * v1 - 2 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3 * f14 * f33 * v2 - 3 * f13 * f34 * v2 +
            f12 * f35 * v2 - f15 * f22 * v3 + 2 * f14 * f23 * v3 - f13 * f24 * v3 +
            2 * f15 * f2 * f3 * v3 - f14 * f22 * f3 * v3 - 4 * f13 * f23 * f3 * v3 +
            3 * f12 * f24 * f3 * v3 - 4 * f14 * f2 * f32 * v3 +
            8 * f13 * f22 * f32 * v3 - 4 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta1 = -(
        (
            -(d2 * f15 * f22) + 2.0 * d2 * f14 * f23 - d2 * f13 * f24 -
            d2 * f14 * f22 * f3 +
            2.0 * d1 * f13 * f23 * f3 +
            2.0 * d2 * f13 * f23 * f3 - 2 * d1 * f12 * f24 * f3 - d2 * f12 * f24 * f3 +
            d2 * f15 * f32 - 3 * d1 * f13 * f22 * f32 - d2 * f13 * f22 * f32 +
            2 * d1 * f12 * f23 * f32 - 2 * d2 * f12 * f23 * f32 +
            d1 * f1 * f24 * f32 +
            2 * d2 * f1 * f24 * f32 - d2 * f14 * f33 +
            d1 * f12 * f22 * f33 +
            3 * d2 * f12 * f22 * f33 - 2 * d1 * f1 * f23 * f33 -
            2 * d2 * f1 * f23 * f33 +
            d1 * f24 * f33 +
            d1 * f13 * f34 +
            d1 * f1 * f22 * f34 - 2 * d1 * f23 * f34 - d1 * f12 * f35 + d1 * f22 * f35 -
            8 * f12 * f23 * f3 * v1 +
            6 * f1 * f24 * f3 * v1 +
            12 * f12 * f22 * f32 * v1 - 8 * f1 * f23 * f32 * v1 - 4 * f12 * f34 * v1 +
            2 * f1 * f35 * v1 +
            2 * f15 * f3 * v2 - 4 * f14 * f32 * v2 + 4 * f12 * f34 * v2 -
            2 * f1 * f35 * v2 - 2 * f15 * f3 * v3 + 8 * f12 * f23 * f3 * v3 -
            6 * f1 * f24 * f3 * v3 + 4 * f14 * f32 * v3 - 12 * f12 * f22 * f32 * v3 +
            8 * f1 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta2 = -(
        (
            d2 * f15 * f2 - d1 * f13 * f23 - 3 * d2 * f13 * f23 +
            d1 * f12 * f24 +
            2.0 * d2 * f12 * f24 - d2 * f15 * f3 + d2 * f14 * f2 * f3 -
            d1 * f12 * f23 * f3 +
            d2 * f12 * f23 * f3 +
            d1 * f1 * f24 * f3 - d2 * f1 * f24 * f3 - d2 * f14 * f32 +
            3 * d1 * f13 * f2 * f32 +
            d2 * f13 * f2 * f32 - d1 * f1 * f23 * f32 + d2 * f1 * f23 * f32 -
            2 * d1 * f24 * f32 - d2 * f24 * f32 - 2 * d1 * f13 * f33 +
            2 * d2 * f13 * f33 - d1 * f12 * f2 * f33 - 3 * d2 * f12 * f2 * f33 +
            3 * d1 * f23 * f33 +
            d2 * f23 * f33 +
            d1 * f12 * f34 - d1 * f1 * f2 * f34 + d1 * f1 * f35 - d1 * f2 * f35 +
            4 * f12 * f23 * v1 - 3 * f1 * f24 * v1 + 4 * f1 * f23 * f3 * v1 -
            3 * f24 * f3 * v1 - 12 * f12 * f2 * f32 * v1 +
            4 * f23 * f32 * v1 +
            8 * f12 * f33 * v1 - f1 * f34 * v1 - f35 * v1 - f15 * v2 - f14 * f3 * v2 +
            8 * f13 * f32 * v2 - 8 * f12 * f33 * v2 +
            f1 * f34 * v2 +
            f35 * v2 +
            f15 * v3 - 4 * f12 * f23 * v3 +
            3 * f1 * f24 * v3 +
            f14 * f3 * v3 - 4 * f1 * f23 * f3 * v3 + 3 * f24 * f3 * v3 -
            8 * f13 * f32 * v3 + 12 * f12 * f2 * f32 * v3 - 4 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta3 = -(
        (
            -2.0 * d2 * f14 * f2 + d1 * f13 * f22 + 3 * d2 * f13 * f22 - d1 * f1 * f24 -
            d2 * f1 * f24 + 2 * d2 * f14 * f3 - 2.0 * d1 * f13 * f2 * f3 -
            2 * d2 * f13 * f2 * f3 + d1 * f12 * f22 * f3 - d2 * f12 * f22 * f3 +
            d1 * f24 * f3 +
            d2 * f24 * f3 +
            d1 * f13 * f32 - d2 * f13 * f32 - 2 * d1 * f12 * f2 * f32 +
            2 * d2 * f12 * f2 * f32 +
            d1 * f1 * f22 * f32 - d2 * f1 * f22 * f32 + d1 * f12 * f33 -
            d2 * f12 * f33 +
            2 * d1 * f1 * f2 * f33 +
            2 * d2 * f1 * f2 * f33 - 3 * d1 * f22 * f33 - d2 * f22 * f33 -
            2 * d1 * f1 * f34 + 2 * d1 * f2 * f34 - 4 * f12 * f22 * v1 +
            2 * f24 * v1 +
            8 * f12 * f2 * f3 * v1 - 4 * f1 * f22 * f3 * v1 - 4 * f12 * f32 * v1 +
            8 * f1 * f2 * f32 * v1 - 4 * f22 * f32 * v1 - 4 * f1 * f33 * v1 +
            2 * f34 * v1 +
            2 * f14 * v2 - 4 * f13 * f3 * v2 + 4 * f1 * f33 * v2 - 2 * f34 * v2 -
            2 * f14 * v3 + 4 * f12 * f22 * v3 - 2 * f24 * v3 + 4 * f13 * f3 * v3 -
            8 * f12 * f2 * f3 * v3 +
            4 * f1 * f22 * f3 * v3 +
            4 * f12 * f32 * v3 - 8 * f1 * f2 * f32 * v3 + 4 * f22 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta4 = -(
        (
            d2 * f13 * f2 - d1 * f12 * f22 - 2 * d2 * f12 * f22 +
            d1 * f1 * f23 +
            d2 * f1 * f23 - d2 * f13 * f3 +
            2.0 * d1 * f12 * f2 * f3 +
            d2 * f12 * f2 * f3 - d1 * f1 * f22 * f3 + d2 * f1 * f22 * f3 -
            d1 * f23 * f3 - d2 * f23 * f3 - d1 * f12 * f32 + d2 * f12 * f32 -
            d1 * f1 * f2 * f32 - 2 * d2 * f1 * f2 * f32 +
            2 * d1 * f22 * f32 +
            d2 * f22 * f32 +
            d1 * f1 * f33 - d1 * f2 * f33 + 3 * f1 * f22 * v1 - 2 * f23 * v1 -
            6 * f1 * f2 * f3 * v1 +
            3 * f22 * f3 * v1 +
            3 * f1 * f32 * v1 - f33 * v1 - f13 * v2 + 3 * f12 * f3 * v2 -
            3 * f1 * f32 * v2 +
            f33 * v2 +
            f13 * v3 - 3 * f1 * f22 * v3 + 2 * f23 * v3 - 3 * f12 * f3 * v3 +
            6 * f1 * f2 * f3 * v3 - 3 * f22 * f3 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )


    # Defined as in LALSimulation - LALSimIMRPhenomUtils.c line 70. Final units are correctly Hz^-1
    # there is a 2 * sqrt(5/(64*pi)) missing w.r.t the standard coefficient, which comes from the (2,2) shperical harmonic
    Overallamp = M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL

    tmpResults = Vector{Any}(undef, 6)
    i = 1
    for ell in (2:4)
        for mm in (ell - 1, ell)
            # Compute ringdown and damping frequencies from fits for the various modes
            fringlm, fdamplm = _RDfreqCalc(model, finMass, aeff, ell, mm)
            Rholm, Taulm = fring / fringlm, fdamplm / fdamp

            # Scale input frequencies according to PhenomHM model
            # Compute mapping coefficinets
            Map_fl = fInsJoin_Ampl
            Map_fi = Map_fl / Rholm
            Map_fr = fringlm

            Map_ai, Map_bi = 2.0 / mm, 0.0

            Map_Trd = Map_fr - fringlm + fring
            Map_Ti = 2.0 * Map_fi / mm
            Map_am = (Map_Trd - Map_Ti) / (Map_fr - Map_fi)
            Map_bm = Map_Ti - Map_fi * Map_am

            Map_ar, Map_br = 1.0, -Map_fr + fring

            # Now scale as f -> f*a+b for each regime
            fgridScaled = @. ifelse(
                fgrid < Map_fi,
                fgrid * Map_ai + Map_bi,
                ifelse(fgrid < Map_fr, fgrid * Map_am + Map_bm, fgrid * Map_ar + Map_br),
            )

            # Map the ampliude's range
            # We divide by the leading order l=m=2 behavior, and then scale in the expected PN behavior for the multipole of interest.

            beta_term1 = _onePointFiveSpinPN_Ampl(fgrid, ell, mm, chi_s, chi_a, eta, Seta)
            beta_term2 = _onePointFiveSpinPN_Ampl(
                2.0 .* fgrid ./ mm,
                ell,
                mm,
                chi_s,
                chi_a,
                eta,
                Seta,
            )

            HMamp_term1 =
                _onePointFiveSpinPN_Ampl(fgridScaled, ell, mm, chi_s, chi_a, eta, Seta)
            HMamp_term2 = _onePointFiveSpinPN_Ampl(fgridScaled, 2, 2, 0.0, 0.0, eta, Seta)

            # The (3,3) and (4,3) modes vanish if eta=0.25 (equal mass case) and the (2,1) mode vanishes if both eta=0.25 and chi1z=chi2z
            completeAmpl = @. Overallamp *
               amp0 *
               (fgridScaled^(-7.0 / 6.0)) *
               ifelse(
                   fgridScaled < fInsJoin_Ampl,
                   1.0 +
                   (fgridScaled^(2.0 / 3.0)) * Acoeffs.two_thirds +
                   (fgridScaled^(4.0 / 3.0)) * Acoeffs.four_thirds +
                   (fgridScaled^(5.0 / 3.0)) * Acoeffs.five_thirds +
                   (fgridScaled^(7.0 / 3.0)) * Acoeffs.seven_thirds +
                   (fgridScaled^(8.0 / 3.0)) * Acoeffs.eight_thirds +
                   fgridScaled * (
                       Acoeffs.one +
                       fgridScaled * Acoeffs.two +
                       fgridScaled * fgridScaled * Acoeffs.three
                   ),
                   ifelse(
                       fgridScaled < fpeak,
                       delta0 +
                       fgridScaled * delta1 +
                       fgridScaled *
                       fgridScaled *
                       (delta2 + fgridScaled * delta3 + fgridScaled * fgridScaled * delta4),
                       ifelse(
                           fgridScaled < fcutPar,
                           exp(-(fgridScaled - fring) * gamma2 / (fdamp * gamma3)) *
                           (fdamp * gamma3 * gamma1) / (
                               (fgridScaled - fring) * (fgridScaled - fring) +
                               fdamp * gamma3 * fdamp * gamma3
                           ),
                           0.0,
                       ),
                   ),
               )

            # This results in NaNs having 0/0, correct for this using nan_to_num()
            tmpResults[i] = replace!(
                completeAmpl .* (beta_term1 ./ beta_term2) .* HMamp_term1 ./ HMamp_term2,
                NaN => 0.0,
            )
            i += 1
        end
    end
    ampllm = AmplslmpStructure(
        tmpResults[1],
        tmpResults[2],
        tmpResults[3],
        tmpResults[4],
        tmpResults[5],
        tmpResults[6],
    )
    return ampllm
end




"""
Compute plus and cross polarisations of the GW as a function of frequency, given the events parameters. The function gives as output h_+ and h_x, which result from
the sum of the contributions of the different modes of the waveform. The function avoids for loops over the modes, dealing with both phase and amplitude in a single function.

    hphc(PhenomHM(), f, mc, eta, chi1, chi2, dL, iota)

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts an array of frequencies.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `dL`: Luminosity distance to the source, in Gpc.
-  `iota`: Inclination angle of the binary. In radians.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
#### Optional arguments:
-  `fInsJoin_PHI`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fInsJoin_Ampl`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate phase. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  (::Tuple{Vector{ComplexF64}, Vector{ComplexF64}}, so a tuple of the two polarisation) GW polarisations h_+ and h_x, as a function of frequency.
#### Example:
```julia
    mc = 30.
    eta = 0.25
    dL = 8.
    iota = pi/3
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomHM(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    hp, hc = hphc(PhenomHM(), f, mc, eta, chi1, chi2, dL, iota)
```
"""

function hphc(model::PhenomHM,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota;
    fcutPar = 0.2,
    fInsJoin_PHI = 0.018,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)

    # This function retuns directly the full plus and cross polarisations, avoiding for loops over the modes

    complShiftm = [0.0, pi * 0.5, 0.0, -pi * 0.5, pi, pi * 0.5, 0.0]

    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    etaInv = 1 ./ eta
    pi2 = pi * pi

    QuadMon1, QuadMon2 = 1.0, 1.0

    chi12, chi22 = chi1 * chi1, chi2 * chi2
    chi1dotchi2 = chi1 * chi2
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    q = 0.5 * (1.0 + Seta - 2.0 * eta) * etaInv
    chi_s2, chi_a2 = chi_s * chi_s, chi_a * chi_a
    chi1dotchi2 = chi1 * chi2
    chi_sdotchi_a = chi_s * chi_a
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # This is MfRef, needed to recover LAL, which sets fRef to f_min if fRef=0
    fRef = fgrid[1] #minimum(fgrid)
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = -1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin, radiated energy and mass
    aeff = _finalspin(model, eta, chi1, chi2)
    Erad = _radiatednrg(model, eta, chi1, chi2)
    finMass = 1.0 - Erad


    # Compute the real and imag parts of the complex ringdown frequency for the (l,m) mode as in LALSimIMRPhenomHM.c line 189
    # These are all fits of the different modes. We directly exploit the fact that the relevant HM in this WF are 6
    modes = [21, 22, 32, 33, 43, 44]
    ells = [2, 2, 3, 3, 4, 4]
    mms = [1, 2, 2, 3, 3, 4]
    # Domain mapping for dimnesionless BH spin
    alphaRDfr = log(2.0 - aeff) / log(3.0)
    # beta = 1./ (2. + l - abs(m))
    betaRDfr = [1.0 / 3.0, 0.5, 1.0 / 3.0, 0.5, 1.0 / 3.0, 0.5]
    kappaRDfr = alphaRDfr .^ betaRDfr
    kappaRDfr2 = kappaRDfr .* kappaRDfr
    kappaRDfr3 = kappaRDfr .* kappaRDfr2
    kappaRDfr4 = kappaRDfr .* kappaRDfr3

    tmpRDfr = @. ifelse(
        modes == 21,
        0.589113 * exp(0.043525 * 1im) +
        0.18896353 * exp(2.289868 * 1im) * kappaRDfr +
        1.15012965 * exp(5.810057 * 1im) * kappaRDfr2 +
        6.04585476 * exp(2.741967 * 1im) * kappaRDfr3 +
        11.12627777 * exp(5.844130 * 1im) * kappaRDfr4 +
        9.34711461 * exp(2.669372 * 1im) * kappaRDfr4 * kappaRDfr +
        3.03838318 * exp(5.791518 * 1im) * kappaRDfr4 * kappaRDfr2,
        ifelse(
            modes == 22,
            1.0 +
            kappaRDfr * (
                1.557847 * exp(2.903124 * 1im) +
                1.95097051 * exp(5.920970 * 1im) * kappaRDfr +
                2.09971716 * exp(2.760585 * 1im) * kappaRDfr2 +
                1.41094660 * exp(5.914340 * 1im) * kappaRDfr3 +
                0.41063923 * exp(2.795235 * 1im) * kappaRDfr4
            ),
            ifelse(
                modes == 32,
                1.022464 * exp(0.004870 * 1im) +
                0.24731213 * exp(0.665292 * 1im) * kappaRDfr +
                1.70468239 * exp(3.138283 * 1im) * kappaRDfr2 +
                0.94604882 * exp(0.163247 * 1im) * kappaRDfr3 +
                1.53189884 * exp(5.703573 * 1im) * kappaRDfr4 +
                2.28052668 * exp(2.685231 * 1im) * kappaRDfr4 * kappaRDfr +
                0.92150314 * exp(5.841704 * 1im) * kappaRDfr4 * kappaRDfr2,
                ifelse(
                    modes == 33,
                    1.5 +
                    kappaRDfr * (
                        2.095657 * exp(2.964973 * 1im) +
                        2.46964352 * exp(5.996734 * 1im) * kappaRDfr +
                        2.66552551 * exp(2.817591 * 1im) * kappaRDfr2 +
                        1.75836443 * exp(5.932693 * 1im) * kappaRDfr3 +
                        0.49905688 * exp(2.781658 * 1im) * kappaRDfr4
                    ),
                    ifelse(
                        modes == 43,
                        1.5 +
                        kappaRDfr * (
                            0.205046 * exp(0.595328 * 1im) +
                            3.10333396 * exp(3.016200 * 1im) * kappaRDfr +
                            4.23612166 * exp(6.038842 * 1im) * kappaRDfr2 +
                            3.02890198 * exp(2.826239 * 1im) * kappaRDfr3 +
                            0.90843949 * exp(5.915164 * 1im) * kappaRDfr4
                        ),
                        2.0 +
                        kappaRDfr * (
                            2.658908 * exp(3.002787 * 1im) +
                            2.97825567 * exp(6.050955 * 1im) * kappaRDfr +
                            3.21842350 * exp(2.877514 * 1im) * kappaRDfr2 +
                            2.12764967 * exp(5.989669 * 1im) * kappaRDfr3 +
                            0.60338186 * exp(2.830031 * 1im) * kappaRDfr4
                        ),
                    ),
                ),
            ),
        ),
    )

    fringlm = real(tmpRDfr) ./ (2.0 * pi * finMass)
    fdamplm = imag(tmpRDfr) ./ (2.0 * pi * finMass)

    # # This recomputation is needed for JAX derivatives
    # betaRDfr = 0.5
    # kappaRDfr  = alphaRDfr^betaRDfr
    # kappaRDfr2 = kappaRDfr*kappaRDfr
    # kappaRDfr3 = kappaRDfr*kappaRDfr2
    # kappaRDfr4 = kappaRDfr*kappaRDfr3

    # tmpRDfr = 1.0 + kappaRDfr * (1.557847 * exp(2.903124 * 1j) + 1.95097051 * exp(5.920970 * 1j) * kappaRDfr + 2.09971716 * exp(2.760585 * 1j) * kappaRDfr2 + 1.41094660 * exp(5.914340 * 1j) * kappaRDfr3 + 0.41063923 * exp(2.795235 * 1j) * kappaRDfr4)

    # fring = (real(tmpRDfr)/(2.*pi*finMass))
    # fdamp = (imag(tmpRDfr)/(2.*pi*finMass))

    #ANDREA: they recalculate the fring and fdamp, but they are already calculated in the previous cell for the 22 mode. I will use the previous ones
    fring = fringlm[2]
    fdamp = fdamplm[2]
    # solved this way
    # Compute sigma coefficients appearing in arXiv:1508.07253 eq. (28)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    sigma1 =
        2096.551999295543 +
        1463.7493168261553 * eta +
        (
            1312.5493286098522 + 18307.330017082117 * eta - 43534.1440746107 * eta2 +
            (-833.2889543511114 + 32047.31997183187 * eta - 108609.45037520859 * eta2) *
            xi +
            (452.25136398112204 + 8353.439546391714 * eta - 44531.3250037322 * eta2) * xi2
        ) * xi
    sigma2 =
        -10114.056472621156 - 44631.01109458185 * eta +
        (
            -6541.308761668722 - 266959.23419307504 * eta +
            686328.3229317984 * eta2 +
            (3405.6372187679685 - 437507.7208209015 * eta + 1.6318171307344697e6 * eta2) *
            xi +
            (-7462.648563007646 - 114585.25177153319 * eta + 674402.4689098676 * eta2) * xi2
        ) * xi
    sigma3 =
        22933.658273436497 +
        230960.00814979506 * eta +
        (
            14961.083974183695 + 1.1940181342318142e6 * eta - 3.1042239693052764e6 * eta2 +
            (-3038.166617199259 + 1.8720322849093592e6 * eta - 7.309145012085539e6 * eta2) *
            xi +
            (42738.22871475411 + 467502.018616601 * eta - 3.064853498512499e6 * eta2) * xi2
        ) * xi
    sigma4 =
        -14621.71522218357 - 377812.8579387104 * eta +
        (
            -9608.682631509726 - 1.7108925257214056e6 * eta +
            4.332924601416521e6 * eta2 +
            (
                -22366.683262266528 - 2.5019716386377467e6 * eta +
                1.0274495902259542e7 * eta2
            ) * xi +
            (-85360.30079034246 - 570025.3441737515 * eta + 4.396844346849777e6 * eta2) *
            xi2
        ) * xi

    # Compute beta coefficients appearing in arXiv:1508.07253 eq. (16)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    beta1 =
        97.89747327985583 - 42.659730877489224 * eta +
        (
            153.48421037904913 - 1417.0620760768954 * eta +
            2752.8614143665027 * eta2 +
            (138.7406469558649 - 1433.6585075135881 * eta + 2857.7418952430758 * eta2) *
            xi +
            (41.025109467376126 - 423.680737974639 * eta + 850.3594335657173 * eta2) * xi2
        ) * xi
    beta2 =
        -3.282701958759534 - 9.051384468245866 * eta +
        (
            -12.415449742258042 + 55.4716447709787 * eta - 106.05109938966335 * eta2 +
            (-11.953044553690658 + 76.80704618365418 * eta - 155.33172948098394 * eta2) *
            xi +
            (-3.4129261592393263 + 25.572377569952536 * eta - 54.408036707740465 * eta2) *
            xi2
        ) * xi
    beta3 =
        -0.000025156429818799565 +
        0.000019750256942201327 * eta +
        (
            -0.000018370671469295915 +
            0.000021886317041311973 * eta +
            0.00008250240316860033 * eta2 +
            (
                7.157371250566708e-6 - 0.000055780000112270685 * eta +
                0.00019142082884072178 * eta2
            ) * xi +
            (
                5.447166261464217e-6 - 0.00003220610095021982 * eta +
                0.00007974016714984341 * eta2
            ) * xi2
        ) * xi

    # Compute alpha coefficients appearing in arXiv:1508.07253 eq. (14)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    alpha1 =
        43.31514709695348 +
        638.6332679188081 * eta +
        (
            -32.85768747216059 + 2415.8938269370315 * eta - 5766.875169379177 * eta2 +
            (-61.85459307173841 + 2953.967762459948 * eta - 8986.29057591497 * eta2) * xi +
            (-21.571435779762044 + 981.2158224673428 * eta - 3239.5664895930286 * eta2) *
            xi2
        ) * xi
    alpha2 =
        -0.07020209449091723 - 0.16269798450687084 * eta +
        (
            -0.1872514685185499 + 1.138313650449945 * eta - 2.8334196304430046 * eta2 +
            (-0.17137955686840617 + 1.7197549338119527 * eta - 4.539717148261272 * eta2) *
            xi +
            (-0.049983437357548705 + 0.6062072055948309 * eta - 1.682769616644546 * eta2) *
            xi2
        ) * xi
    alpha3 =
        9.5988072383479 - 397.05438595557433 * eta +
        (
            16.202126189517813 - 1574.8286986717037 * eta +
            3600.3410843831093 * eta2 +
            (27.092429659075467 - 1786.482357315139 * eta + 5152.919378666511 * eta2) * xi +
            (11.175710130033895 - 577.7999423177481 * eta + 1808.730762932043 * eta2) * xi2
        ) * xi
    alpha4 =
        -0.02989487384493607 +
        1.4022106448583738 * eta +
        (
            -0.07356049468633846 +
            0.8337006542278661 * eta +
            0.2240008282397391 * eta2 +
            (-0.055202870001177226 + 0.5667186343606578 * eta + 0.7186931973380503 * eta2) *
            xi +
            (
                -0.015507437354325743 +
                0.15750322779277187 * eta +
                0.21076815715176228 * eta2
            ) * xi2
        ) * xi
    alpha5 =
        0.9974408278363099 - 0.007884449714907203 * eta +
        (
            -0.059046901195591035 + 1.3958712396764088 * eta - 4.516631601676276 * eta2 +
            (-0.05585343136869692 + 1.7516580039343603 * eta - 5.990208965347804 * eta2) *
            xi +
            (-0.017945336522161195 + 0.5965097794825992 * eta - 2.0608879367971804 * eta2) *
            xi2
        ) * xi

    # Compute the TF2 phase coefficients and put them in a dictionary (spin effects are included up to 3.5PN)

    TF2OverallAmpl = 3 / (128.0 * eta)
    TF2_5coeff_tmp =
        38645.0 * pi / 756.0 - 65.0 * pi * eta / 9.0 - (
            (732985.0 / 2268.0 - 24260.0 * eta / 81.0 - 340.0 * eta2 / 9.0) * chi_s +
            (732985.0 / 2268.0 + 140.0 * eta / 9.0) * Seta * chi_a
        ) #variable to be used later
    TF2_6coeff_tmp =
        11583.231236531 / 4.694215680 - 640.0 / 3.0 * pi2 -
        684.8 / 2.1 * MathConstants.eulergamma +
        eta * (-15737.765635 / 3.048192 + 225.5 / 1.2 * pi2) +
        eta2 * 76.055 / 1.728 - eta2 * eta * 127.825 / 1.296 - log(4.0) * 684.8 / 2.1 +
        pi * chi1 * m1ByM * (1490.0 / 3.0 + m1ByM * 260.0) +
        pi * chi2 * m2ByM * (1490.0 / 3.0 + m2ByM * 260.0) +
        (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) *
        m1ByM^2 *
        QuadMon1 *
        chi12 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 ) *
        m1ByM^2 *
        chi12 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) *
        m2ByM^2 *
        QuadMon2 *
        chi22 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 ) *
        m2ByM^2 *
        chi22

    TF2coeffs = TF2coeffsStructure(
        1.0,
        0.0,
        3715.0 / 756.0 + (55.0 * eta) / 9.0,
        -16.0 * pi +
        (113.0 * Seta * chi_a) / 3.0 +
        (113.0 / 3.0 - (76.0 * eta) / 3.0) * chi_s,
        5.0 * (3058.673 / 7.056 + 5429.0 / 7.0 * eta + 617.0 * eta2) / 72.0 +
        247.0 / 4.8 * eta * chi1dotchi2 - 721.0 / 4.8 * eta * chi1dotchi2 +
        (-720.0 / 9.6 * QuadMon1 + 1.0 / 9.6) * m1ByM^2  * chi12 +
        (-720.0 / 9.6 * QuadMon2 + 1.0 / 9.6) * m2ByM^2  * chi22 +
        (240.0 / 9.6 * QuadMon1 - 7.0 / 9.6) * m1ByM^2  * chi12 +
        (240.0 / 9.6 * QuadMon2 - 7.0 / 9.6) * m2ByM^2  * chi22,
        TF2_5coeff_tmp,
        TF2_5coeff_tmp * 3.0,
        TF2_6coeff_tmp - (
            (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 )
            ) *
            m1ByM *
            m1ByM *
            chi12 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 )
            ) *
            m2ByM *
            m2ByM *
            chi22
        ),
        -6848.0 / 21.0,
        77096675.0 * pi / 254016.0 + 378515.0 * pi * eta / 1512.0 -
        74045.0 * pi * eta2 / 756.0 +
        (
            -25150083775.0 / 3048192.0 + 10566655595.0 * eta / 762048.0 -
            1042165.0 * eta2 / 3024.0 + 5345.0 * eta2 * eta / 36.0
        ) * chi_s +
        Seta * (
            (
                -25150083775.0 / 3048192.0 + 26804935.0 * eta / 6048.0 -
                1985.0 * eta2 / 48.0
            ) * chi_a
        ),
    )


    PhiInspcoeffs = PhiInspcoeffsStructure(
        TF2coeffs.five * TF2OverallAmpl,
        TF2coeffs.seven * TF2OverallAmpl * (pi^(2.0 / 3.0)),
        TF2coeffs.six * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.six_log * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.five_log * TF2OverallAmpl,
        TF2coeffs.four * TF2OverallAmpl * (pi^(-1.0 / 3.0)),
        TF2coeffs.three * TF2OverallAmpl * (pi^(-2.0 / 3.0)),
        TF2coeffs.two * TF2OverallAmpl / pi,
        TF2coeffs.one * TF2OverallAmpl * (pi^(-4.0 / 3.0)),
        TF2coeffs.zero * TF2OverallAmpl * (pi^(-5.0 / 3.0)),
        sigma1,
        sigma2 * 0.75,
        sigma3 * 0.6,
        sigma4 * 0.5,
    )
    #Now compute the coefficients to align the three parts

    fInsJoinPh = fInsJoin_PHI
    fMRDJoinPh = 0.5 * fring

    # First the Inspiral - Intermediate: we compute C1Int and C2Int coeffs
    # Equations to solve for to get C(1) continuous join
    # PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    # Joining at fInsJoin
    # PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    # PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    # This is the first derivative wrt f of the inspiral phase computed at fInsJoin, first add the PN contribution and then the higher order calibrated terms
    DPhiIns =
        (
            2.0 * TF2coeffs.seven * TF2OverallAmpl * ((pi * fInsJoinPh)^(7.0 / 3.0)) +
            (
                TF2coeffs.six * TF2OverallAmpl +
                TF2coeffs.six_log * TF2OverallAmpl * (1.0 + log(pi * fInsJoinPh) / 3.0)
            ) * ((pi * fInsJoinPh)^(2.0)) +
            TF2coeffs.five_log * TF2OverallAmpl * ((pi * fInsJoinPh)^(5.0 / 3.0)) -
            TF2coeffs.four * TF2OverallAmpl * ((pi * fInsJoinPh)^(4.0 / 3.0)) -
            2.0 * TF2coeffs.three * TF2OverallAmpl * (pi * fInsJoinPh) -
            3.0 * TF2coeffs.two * TF2OverallAmpl * ((pi * fInsJoinPh)^(2.0 / 3.0)) -
            4.0 * TF2coeffs.one * TF2OverallAmpl * ((pi * fInsJoinPh)^(1.0 / 3.0)) -
            5.0 * TF2coeffs.zero * TF2OverallAmpl
        ) * pi / (3.0 * ((pi * fInsJoinPh)^(8.0 / 3.0))) +
        (
            sigma1 +
            sigma2 * (fInsJoinPh^(1.0 / 3.0)) +
            sigma3 * (fInsJoinPh^(2.0 / 3.0)) +
            sigma4 * fInsJoinPh
        ) * etaInv
    # This is the first derivative of the Intermediate phase computed at fMRDJoin
    DPhiInt = (beta1 + beta3 / (fInsJoinPh^4) + beta2 / fInsJoinPh) * etaInv
    C2Int = DPhiIns - DPhiInt

    # This is the inspiral phase computed at fInsJoin
    PhiInsJoin =
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (fInsJoinPh^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (fInsJoinPh^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (fInsJoinPh^(1.0 / 3.0)) * log(pi * fInsJoinPh) / 3.0 +
        PhiInspcoeffs.log * log(pi * fInsJoinPh) / 3.0 +
        PhiInspcoeffs.min_third * (fInsJoinPh^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (fInsJoinPh^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / fInsJoinPh +
        PhiInspcoeffs.min_four_thirds * (fInsJoinPh^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (fInsJoinPh^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * fInsJoinPh +
            PhiInspcoeffs.four_thirds * (fInsJoinPh^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (fInsJoinPh^(5.0 / 3.0)) +
            PhiInspcoeffs.two * fInsJoinPh * fInsJoinPh
        ) * etaInv
    # This is the Intermediate phase computed at fMRDJoinPh
    PhiIntJoin =
        beta1 * fInsJoinPh - beta3 / (3.0 * fInsJoinPh * fInsJoinPh * fInsJoinPh) +
        beta2 * log(fInsJoinPh)

    C1Int = PhiInsJoin - PhiIntJoin * etaInv - C2Int * fInsJoinPh
    # Now the same for Intermediate - Merger-Ringdown: we also need a temporary Intermediate Phase function
    PhiIntTempVal =
        (beta1 * fMRDJoinPh - beta3 / (3.0 * fMRDJoinPh^3) + beta2 * log(fMRDJoinPh)) *
        etaInv +
        C1Int +
        C2Int * fMRDJoinPh
    DPhiIntTempVal = C2Int + (beta1 + beta3 / (fMRDJoinPh^4) + beta2 / fMRDJoinPh) * etaInv

    DPhiMRDVal =
        (
            alpha1 +
            alpha2 / (fMRDJoinPh^2) +
            alpha3 / (fMRDJoinPh^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fMRDJoinPh - alpha5 * fring) * (fMRDJoinPh - alpha5 * fring) /
                    (fdamp^2)
                )
            )
        ) * etaInv
    PhiMRJoinTemp =
        -(alpha2 / fMRDJoinPh) +
        (4.0 / 3.0) * (alpha3 * (fMRDJoinPh^(0.75))) +
        alpha1 * fMRDJoinPh +
        alpha4 * atan((fMRDJoinPh - alpha5 * fring) / fdamp)
    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRD * fMRDJoinPh

    # Compute coefficients gamma appearing in arXiv:1508.07253 eq. (19), the numerical coefficients are in Tab. 5
    gamma1 =
        0.006927402739328343 +
        0.03020474290328911 * eta +
        (
            0.006308024337706171 - 0.12074130661131138 * eta +
            0.26271598905781324 * eta2 +
            (
                0.0034151773647198794 - 0.10779338611188374 * eta +
                0.27098966966891747 * eta2
            ) * xi +
            (
                0.0007374185938559283 - 0.02749621038376281 * eta +
                0.0733150789135702 * eta2
            ) * xi2
        ) * xi
    gamma2 =
        1.010344404799477 +
        0.0008993122007234548 * eta +
        (
            0.283949116804459 - 4.049752962958005 * eta +
            13.207828172665366 * eta2 +
            (0.10396278486805426 - 7.025059158961947 * eta + 24.784892370130475 * eta2) *
            xi +
            (0.03093202475605892 - 2.6924023896851663 * eta + 9.609374464684983 * eta2) *
            xi2
        ) * xi
    gamma3 =
        1.3081615607036106 - 0.005537729694807678 * eta +
        (
            -0.06782917938621007 - 0.6689834970767117 * eta +
            3.403147966134083 * eta2 +
            (-0.05296577374411866 - 0.9923793203111362 * eta + 4.820681208409587 * eta2) *
            xi +
            (
                -0.006134139870393713 - 0.38429253308696365 * eta +
                1.7561754421985984 * eta2
            ) * xi2
        ) * xi
    # Compute fpeak, from arXiv:1508.07253 eq. (20), we remove the square root term in case it is complex
    fpeak = ifelse(
        gamma2 >= 1.0,
        abs(fring - (fdamp * gamma3) / gamma2),
        fring + (fdamp * (-1.0 + sqrt(1.0 - gamma2 * gamma2)) * gamma3) / gamma2,
    )
    # Compute coefficients rho appearing in arXiv:1508.07253 eq. (30), the numerical coefficients are in Tab. 5
    rho1 =
        3931.8979897196696 - 17395.758706812805 * eta +
        (
            3132.375545898835 + 343965.86092361377 * eta - 1.2162565819981997e6 * eta2 +
            (-70698.00600428853 + 1.383907177859705e6 * eta - 3.9662761890979446e6 * eta2) *
            xi +
            (-60017.52423652596 + 803515.1181825735 * eta - 2.091710365941658e6 * eta2) *
            xi2
        ) * xi
    rho2 =
        -40105.47653771657 +
        112253.0169706701 * eta +
        (
            23561.696065836168 - 3.476180699403351e6 * eta +
            1.137593670849482e7 * eta2 +
            (754313.1127166454 - 1.308476044625268e7 * eta + 3.6444584853928134e7 * eta2) *
            xi +
            (596226.612472288 - 7.4277901143564405e6 * eta + 1.8928977514040343e7 * eta2) *
            xi2
        ) * xi
    rho3 =
        83208.35471266537 - 191237.7264145924 * eta +
        (
            -210916.2454782992 + 8.71797508352568e6 * eta - 2.6914942420669552e7 * eta2 +
            (
                -1.9889806527362722e6 + 3.0888029960154563e7 * eta -
                8.390870279256162e7 * eta2
            ) * xi +
            (
                -1.4535031953446497e6 + 1.7063528990822166e7 * eta -
                4.2748659731120914e7 * eta2
            ) * xi2
        ) * xi
    # Compute coefficients delta appearing in arXiv:1508.07253 eq. (21)
    f1Interm = fInsJoin_Ampl
    f3Interm = fpeak
    dfInterm = 0.5 * (f3Interm - f1Interm)
    f2Interm = f1Interm + dfInterm
    # First write the inspiral coefficients, we put them in a dictionary and label with the power in front of which they appear
    amp0 = sqrt(2.0 * eta / 3.0) * (pi^(-1.0 / 6.0))
    Acoeffs = AcoeffsStructure(
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / 672.0,
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0,
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) * (pi^(4.0 / 3.0))
        ) / 8.128512e6,
        (
            (pi^(5.0 / 3.0)) * (
                chi2 * (
                    -285197.0 * (-1.0 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1.0 + 4.0 * eta) * pi
            )
        ) / 32256.0,
        -(
            (pi^2.0) * (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368.0 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320.0 * pi2
                )
            )
        ) / 6.0085960704e10,
        rho1,
        rho2,
        rho3,
    )

    # v1 is the inspiral model evaluated at f1Interm
    v1 =
        1.0 +
        (f1Interm^(2.0 / 3.0)) * Acoeffs.two_thirds +
        (f1Interm^(4.0 / 3.0)) * Acoeffs.four_thirds +
        (f1Interm^(5.0 / 3.0)) * Acoeffs.five_thirds +
        (f1Interm^(7.0 / 3.0)) * Acoeffs.seven_thirds +
        (f1Interm^(8.0 / 3.0)) * Acoeffs.eight_thirds +
        f1Interm *
        (Acoeffs.one + f1Interm * Acoeffs.two + f1Interm * f1Interm * Acoeffs.three)
    # d1 is the derivative of the inspiral model evaluated at f1
    d1 =
        ((-969.0 + 1804.0 * eta) * (pi^(2.0 / 3.0))) / (1008.0 * (f1Interm^(1.0 / 3.0))) +
        (
            (
                chi1 * (81.0 * SetaPlus1 - 44.0 * eta) +
                chi2 * (81.0 - 81.0 * Seta - 44.0 * eta)
            ) * pi
        ) / 48.0 +
        (
            (
                -27312085.0 - 10287648.0 * chi22 - 10287648.0 * chi12 * SetaPlus1 +
                10287648.0 * chi22 * Seta +
                24.0 *
                (
                    -1975055.0 + 857304.0 * chi12 - 994896.0 * chi1 * chi2 +
                    857304.0 * chi22
                ) *
                eta +
                35371056.0 * eta2
            ) *
            (f1Interm^(1.0 / 3.0)) *
            (pi^(4.0 / 3.0))
        ) / 6.096384e6 +
        (
            5.0 *
            (f1Interm^(2.0 / 3.0)) *
            (pi^(5.0 / 3.0)) *
            (
                chi2 * (
                    -285197.0 * (-1 + Seta) + 4.0 * (-91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                chi1 * (
                    285197.0 * SetaPlus1 - 4.0 * (91902.0 + 1579.0 * Seta) * eta -
                    35632.0 * eta2
                ) +
                42840.0 * (-1 + 4 * eta) * pi
            )
        ) / 96768.0 -
        (
            f1Interm *
            pi2 *
            (
                -336.0 *
                (
                    -3248849057.0 + 2943675504.0 * chi12 - 3339284256.0 * chi1 * chi2 +
                    2943675504.0 * chi22
                ) *
                eta2 - 324322727232.0 * eta2 * eta -
                7.0 * (
                    -177520268561.0 +
                    107414046432.0 * chi22 +
                    107414046432.0 * chi12 * SetaPlus1 - 107414046432.0 * chi22 * Seta +
                    11087290368 * (chi1 + chi2 + chi1 * Seta - chi2 * Seta) * pi
                ) +
                12.0 *
                eta *
                (
                    -545384828789.0 - 176491177632.0 * chi1 * chi2 +
                    202603761360.0 * chi22 +
                    77616.0 * chi12 * (2610335.0 + 995766.0 * Seta) -
                    77287373856.0 * chi22 * Seta +
                    5841690624.0 * (chi1 + chi2) * pi +
                    21384760320 * pi2
                )
            )
        ) / 3.0042980352e10 +
        (7.0 / 3.0) * (f1Interm^(4.0 / 3.0)) * rho1 +
        (8.0 / 3.0) * (f1Interm^(5.0 / 3.0)) * rho2 +
        3.0 * (f1Interm * f1Interm) * rho3
    # v3 is the merger-ringdown model (eq. (19) of arXiv:1508.07253) evaluated at f3
    v3 =
        exp(-(f3Interm - fring) * gamma2 / (fdamp * gamma3)) * (fdamp * gamma3 * gamma1) /
        ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
    # d2 is the derivative of the merger-ringdown model evaluated at f3
    d2 =
        (
            (-2.0 * fdamp * (f3Interm - fring) * gamma3 * gamma1) /
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3) -
            (gamma2 * gamma1)
        ) / (
            exp((f3Interm - fring) * gamma2 / (fdamp * gamma3)) *
            ((f3Interm - fring) * (f3Interm - fring) + fdamp * gamma3 * fdamp * gamma3)
        )
    # v2 is the value of the amplitude evaluated at f2. They come from the fit of the collocation points in the intermediate region
    v2 =
        0.8149838730507785 +
        2.5747553517454658 * eta +
        (
            1.1610198035496786 - 2.3627771785551537 * eta +
            6.771038707057573 * eta2 +
            (0.7570782938606834 - 2.7256896890432474 * eta + 7.1140380397149965 * eta2) *
            xi +
            (0.1766934149293479 - 0.7978690983168183 * eta + 2.1162391502005153 * eta2) *
            xi2
        ) * xi
    # Now some definitions to speed up
    f1 = f1Interm
    f2 = f2Interm
    f3 = f3Interm
    f12 = f1Interm * f1Interm
    f13 = f1Interm * f12
    f14 = f1Interm * f13
    f15 = f1Interm * f14
    f22 = f2Interm * f2Interm
    f23 = f2Interm * f22
    f24 = f2Interm * f23
    f32 = f3Interm * f3Interm
    f33 = f3Interm * f32
    f34 = f3Interm * f33
    f35 = f3Interm * f34
    # Finally conpute the deltas
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2.0 * d1 * f13 * f22 * f33 - 2.0 * d2 * f13 * f22 * f33 -
            d1 * f12 * f23 * f33 + d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 -
            d1 * f13 * f2 * f34 - d1 * f12 * f22 * f34 +
            2.0 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4.0 * f12 * f23 * f32 * v1 -
            3.0 * f1 * f24 * f32 * v1 - 8.0 * f12 * f22 * f33 * v1 +
            4.0 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4.0 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2.0 * f23 * f34 * v1 - 2.0 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3.0 * f14 * f33 * v2 -
            3.0 * f13 * f34 * v2 + f12 * f35 * v2 - f15 * f22 * v3 +
            2.0 * f14 * f23 * v3 - f13 * f24 * v3 + 2.0 * f15 * f2 * f3 * v3 -
            f14 * f22 * f3 * v3 - 4.0 * f13 * f23 * f3 * v3 +
            3.0 * f12 * f24 * f3 * v3 - 4.0 * f14 * f2 * f32 * v3 +
            8.0 * f13 * f22 * f32 * v3 - 4.0 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta0 = -(
        (
            d2 * f15 * f22 * f3 - 2.0 * d2 * f14 * f23 * f3 + d2 * f13 * f24 * f3 -
            d2 * f15 * f2 * f32 + d2 * f14 * f22 * f32 - d1 * f13 * f23 * f32 +
            d2 * f13 * f23 * f32 +
            d1 * f12 * f24 * f32 - d2 * f12 * f24 * f32 +
            d2 * f14 * f2 * f33 +
            2 * d1 * f13 * f22 * f33 - 2 * d2 * f13 * f22 * f33 - d1 * f12 * f23 * f33 +
            d2 * f12 * f23 * f33 - d1 * f1 * f24 * f33 - d1 * f13 * f2 * f34 -
            d1 * f12 * f22 * f34 +
            2 * d1 * f1 * f23 * f34 +
            d1 * f12 * f2 * f35 - d1 * f1 * f22 * f35 + 4 * f12 * f23 * f32 * v1 -
            3 * f1 * f24 * f32 * v1 - 8 * f12 * f22 * f33 * v1 +
            4 * f1 * f23 * f33 * v1 +
            f24 * f33 * v1 +
            4 * f12 * f2 * f34 * v1 +
            f1 * f22 * f34 * v1 - 2 * f23 * f34 * v1 - 2 * f1 * f2 * f35 * v1 +
            f22 * f35 * v1 - f15 * f32 * v2 + 3 * f14 * f33 * v2 - 3 * f13 * f34 * v2 +
            f12 * f35 * v2 - f15 * f22 * v3 + 2 * f14 * f23 * v3 - f13 * f24 * v3 +
            2 * f15 * f2 * f3 * v3 - f14 * f22 * f3 * v3 - 4 * f13 * f23 * f3 * v3 +
            3 * f12 * f24 * f3 * v3 - 4 * f14 * f2 * f32 * v3 +
            8 * f13 * f22 * f32 * v3 - 4 * f12 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta1 = -(
        (
            -(d2 * f15 * f22) + 2.0 * d2 * f14 * f23 - d2 * f13 * f24 -
            d2 * f14 * f22 * f3 +
            2.0 * d1 * f13 * f23 * f3 +
            2.0 * d2 * f13 * f23 * f3 - 2 * d1 * f12 * f24 * f3 - d2 * f12 * f24 * f3 +
            d2 * f15 * f32 - 3 * d1 * f13 * f22 * f32 - d2 * f13 * f22 * f32 +
            2 * d1 * f12 * f23 * f32 - 2 * d2 * f12 * f23 * f32 +
            d1 * f1 * f24 * f32 +
            2 * d2 * f1 * f24 * f32 - d2 * f14 * f33 +
            d1 * f12 * f22 * f33 +
            3 * d2 * f12 * f22 * f33 - 2 * d1 * f1 * f23 * f33 -
            2 * d2 * f1 * f23 * f33 +
            d1 * f24 * f33 +
            d1 * f13 * f34 +
            d1 * f1 * f22 * f34 - 2 * d1 * f23 * f34 - d1 * f12 * f35 + d1 * f22 * f35 -
            8 * f12 * f23 * f3 * v1 +
            6 * f1 * f24 * f3 * v1 +
            12 * f12 * f22 * f32 * v1 - 8 * f1 * f23 * f32 * v1 - 4 * f12 * f34 * v1 +
            2 * f1 * f35 * v1 +
            2 * f15 * f3 * v2 - 4 * f14 * f32 * v2 + 4 * f12 * f34 * v2 -
            2 * f1 * f35 * v2 - 2 * f15 * f3 * v3 + 8 * f12 * f23 * f3 * v3 -
            6 * f1 * f24 * f3 * v3 + 4 * f14 * f32 * v3 - 12 * f12 * f22 * f32 * v3 +
            8 * f1 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta2 = -(
        (
            d2 * f15 * f2 - d1 * f13 * f23 - 3 * d2 * f13 * f23 +
            d1 * f12 * f24 +
            2.0 * d2 * f12 * f24 - d2 * f15 * f3 + d2 * f14 * f2 * f3 -
            d1 * f12 * f23 * f3 +
            d2 * f12 * f23 * f3 +
            d1 * f1 * f24 * f3 - d2 * f1 * f24 * f3 - d2 * f14 * f32 +
            3 * d1 * f13 * f2 * f32 +
            d2 * f13 * f2 * f32 - d1 * f1 * f23 * f32 + d2 * f1 * f23 * f32 -
            2 * d1 * f24 * f32 - d2 * f24 * f32 - 2 * d1 * f13 * f33 +
            2 * d2 * f13 * f33 - d1 * f12 * f2 * f33 - 3 * d2 * f12 * f2 * f33 +
            3 * d1 * f23 * f33 +
            d2 * f23 * f33 +
            d1 * f12 * f34 - d1 * f1 * f2 * f34 + d1 * f1 * f35 - d1 * f2 * f35 +
            4 * f12 * f23 * v1 - 3 * f1 * f24 * v1 + 4 * f1 * f23 * f3 * v1 -
            3 * f24 * f3 * v1 - 12 * f12 * f2 * f32 * v1 +
            4 * f23 * f32 * v1 +
            8 * f12 * f33 * v1 - f1 * f34 * v1 - f35 * v1 - f15 * v2 - f14 * f3 * v2 +
            8 * f13 * f32 * v2 - 8 * f12 * f33 * v2 +
            f1 * f34 * v2 +
            f35 * v2 +
            f15 * v3 - 4 * f12 * f23 * v3 +
            3 * f1 * f24 * v3 +
            f14 * f3 * v3 - 4 * f1 * f23 * f3 * v3 + 3 * f24 * f3 * v3 -
            8 * f13 * f32 * v3 + 12 * f12 * f2 * f32 * v3 - 4 * f23 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta3 = -(
        (
            -2.0 * d2 * f14 * f2 + d1 * f13 * f22 + 3 * d2 * f13 * f22 - d1 * f1 * f24 -
            d2 * f1 * f24 + 2 * d2 * f14 * f3 - 2.0 * d1 * f13 * f2 * f3 -
            2 * d2 * f13 * f2 * f3 + d1 * f12 * f22 * f3 - d2 * f12 * f22 * f3 +
            d1 * f24 * f3 +
            d2 * f24 * f3 +
            d1 * f13 * f32 - d2 * f13 * f32 - 2 * d1 * f12 * f2 * f32 +
            2 * d2 * f12 * f2 * f32 +
            d1 * f1 * f22 * f32 - d2 * f1 * f22 * f32 + d1 * f12 * f33 -
            d2 * f12 * f33 +
            2 * d1 * f1 * f2 * f33 +
            2 * d2 * f1 * f2 * f33 - 3 * d1 * f22 * f33 - d2 * f22 * f33 -
            2 * d1 * f1 * f34 + 2 * d1 * f2 * f34 - 4 * f12 * f22 * v1 +
            2 * f24 * v1 +
            8 * f12 * f2 * f3 * v1 - 4 * f1 * f22 * f3 * v1 - 4 * f12 * f32 * v1 +
            8 * f1 * f2 * f32 * v1 - 4 * f22 * f32 * v1 - 4 * f1 * f33 * v1 +
            2 * f34 * v1 +
            2 * f14 * v2 - 4 * f13 * f3 * v2 + 4 * f1 * f33 * v2 - 2 * f34 * v2 -
            2 * f14 * v3 + 4 * f12 * f22 * v3 - 2 * f24 * v3 + 4 * f13 * f3 * v3 -
            8 * f12 * f2 * f3 * v3 +
            4 * f1 * f22 * f3 * v3 +
            4 * f12 * f32 * v3 - 8 * f1 * f2 * f32 * v3 + 4 * f22 * f32 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )
    delta4 = -(
        (
            d2 * f13 * f2 - d1 * f12 * f22 - 2 * d2 * f12 * f22 +
            d1 * f1 * f23 +
            d2 * f1 * f23 - d2 * f13 * f3 +
            2.0 * d1 * f12 * f2 * f3 +
            d2 * f12 * f2 * f3 - d1 * f1 * f22 * f3 + d2 * f1 * f22 * f3 -
            d1 * f23 * f3 - d2 * f23 * f3 - d1 * f12 * f32 + d2 * f12 * f32 -
            d1 * f1 * f2 * f32 - 2 * d2 * f1 * f2 * f32 +
            2 * d1 * f22 * f32 +
            d2 * f22 * f32 +
            d1 * f1 * f33 - d1 * f2 * f33 + 3 * f1 * f22 * v1 - 2 * f23 * v1 -
            6 * f1 * f2 * f3 * v1 +
            3 * f22 * f3 * v1 +
            3 * f1 * f32 * v1 - f33 * v1 - f13 * v2 + 3 * f12 * f3 * v2 -
            3 * f1 * f32 * v2 +
            f33 * v2 +
            f13 * v3 - 3 * f1 * f22 * v3 + 2 * f23 * v3 - 3 * f12 * f3 * v3 +
            6 * f1 * f2 * f3 * v3 - 3 * f22 * f3 * v3
        ) / (
            (f1 - f2) *
            (f1 - f2) *
            (f1 - f3) *
            (f1 - f3) *
            (f1 - f3) *
            (f3 - f2) *
            (f3 - f2)
        )
    )

    # Defined as in LALSimulation - LALSimIMRPhenomUtils.c line 70. Final units are correctly Hz^-1
    # there is a 2 * sqrt(5/(64*pi)) missing w.r.t the standard coefficient, which comes from the (2,2) shperical harmonic
    Overallamp = M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL

    # function completeAmpl(infreqs, Acoeffs=Acoeffs, delta0=delta0, delta1=delta1, delta2=delta2, delta3=delta3, delta4=delta4, fpeak=fpeak, fring=fring, fdamp=fdamp, fInsJoin_Ampl=fInsJoin_Ampl, fcutPar=fcutPar)

    #     return @. Overallamp*amp0*(infreqs^(-7. /6.))*ifelse(infreqs < fInsJoin_Ampl, 1. + (infreqs^(2. /3.))*Acoeffs.two_thirds + (infreqs^(4. /3.))*Acoeffs.four_thirds + (infreqs^(5. /3.))*Acoeffs.five_thirds + (infreqs^(7. /3.))*Acoeffs.seven_thirds + (infreqs^(8. /3.))*Acoeffs.eight_thirds + infreqs*(Acoeffs.one + infreqs*Acoeffs.two + infreqs*infreqs*Acoeffs.three), ifelse(infreqs < fpeak, delta0 + infreqs*delta1 + infreqs*infreqs*(delta2 + infreqs*delta3 + infreqs*infreqs*delta4), ifelse(infreqs < fcutPar, exp(-(infreqs - fring)*gamma2/(fdamp*gamma3))* (fdamp*gamma3*gamma1) / ((infreqs - fring)*(infreqs - fring) + fdamp*gamma3*fdamp*gamma3), 0.)))
    # end




    # Time shift so that peak amplitude is approximately at t=0
    t0 =
        (
            alpha1 +
            alpha2 / (fpeak * fpeak) +
            alpha3 / (fpeak^(1.0 / 4.0)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fpeak - alpha5 * fring) * (fpeak - alpha5 * fring) / (fdamp * fdamp)
                )
            )
        ) * etaInv

    phiRef = _completePhase(
        fRef,
        C1MRD,
        C2MRD,
        1.0,
        1.0,
        C1Int,
        C2Int,
        PhiInspcoeffs,
        fInsJoin_PHI,
        alpha1,
        alpha2,
        alpha3,
        alpha4,
        alpha5,
        beta1,
        beta2,
        beta3,
        etaInv,
        fMRDJoinPh,
        fcutPar,
        fring,
        fdamp,
    )

    phi0 = 0.5 * phiRef

    # Now compute all the modes, they are 6, we parallelize

    Rholm, Taulm = fring ./ fringlm, fdamplm ./ fdamp
    # Rholm and Taulm only figure in the MRD part, the rest of the coefficients is the same, recompute only this
    DPhiMRDVal = @. (
        alpha1 +
        alpha2 / (fMRDJoinPh * fMRDJoinPh) +
        alpha3 / (fMRDJoinPh^(1.0 / 4.0)) +
        alpha4 / (
            fdamp *
            Taulm *
            (
                1.0 +
                (fMRDJoinPh - alpha5 * fring) * (fMRDJoinPh - alpha5 * fring) /
                (fdamp * Taulm * Rholm * fdamp * Taulm * Rholm)
            )
        )
    ) * etaInv
    PhiMRJoinTemp = @. -(alpha2 / fMRDJoinPh) +
       (4.0 / 3.0) * (alpha3 * (fMRDJoinPh^(3.0 / 4.0))) +
       alpha1 * fMRDJoinPh +
       alpha4 * Rholm * atan((fMRDJoinPh - alpha5 * fring) / (fdamp * Rholm * Taulm))
    C2MRDHM = @. DPhiIntTempVal - DPhiMRDVal
    C1MRDHM = @. (PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRDHM * fMRDJoinPh)
    #Rholm, Taulm, DPhiMRDVal, PhiMRJoinTemp, C2MRDHM = Rholm', Taulm', DPhiMRDVal', PhiMRJoinTemp', C2MRDHM'

    # Scale input frequencies according to PhenomHM model
    # Compute mapping coefficinets
    Map_flPhi = fInsJoin_PHI
    Map_fiPhi = Map_flPhi ./ Rholm
    Map_flAmp = fInsJoin_Ampl
    Map_fiAmp = Map_flAmp ./ Rholm
    Map_fr = fringlm

    Map_ai, Map_bi = 2.0 ./ mms, 0.0


    Map_TrdAmp = @. Map_fr - fringlm + fring
    Map_TiAmp = @. 2.0 * Map_fiAmp / mms
    Map_amAmp = @. (Map_TrdAmp - Map_TiAmp) / (Map_fr - Map_fiAmp)
    Map_bmAmp = @. Map_TiAmp - Map_fiAmp * Map_amAmp

    Map_TrdPhi = @. Map_fr * Rholm
    Map_TiPhi = @. 2.0 * Map_fiPhi / mms
    Map_amPhi = @. (Map_TrdPhi - Map_TiPhi) / (Map_fr - Map_fiPhi)
    Map_bmPhi = @. Map_TiPhi - Map_fiPhi * Map_amPhi

    Map_arAmp, Map_brAmp = 1.0, -Map_fr .+ fring
    Map_arPhi, Map_brPhi = Rholm, 0.0

    # Here we do a reshape of the input frequencies, to make the operations easier
    # We convert them to matrices (1-dim matrices) because we need to broadcast the operations
    # For example a 1-dim matrix of 1000x1 elements and a 1-dim matrix of 1x6 elements, when multiplied, will give a 1000x6 matrix
    # It can be hard to follow but remember that the matrices (1000,1) are frequencies or quantities independent of the modes, 
    # and the matrices (1,6) are the coefficients of the modes
    # The result is a matrix (1000,6) with the values of the modes for each frequency
    fgrid = reshape(fgrid, :, 1)
    Map_ai, Map_amAmp, Map_bmAmp, Map_brAmp = reshape(Map_ai, 1, :),
    reshape(Map_amAmp, 1, :),
    reshape(Map_bmAmp, 1, :),
    reshape(Map_brAmp, 1, :)
    Map_amPhi, Map_bmPhi, Map_arPhi =
        reshape(Map_amPhi, 1, :), reshape(Map_bmPhi, 1, :), reshape(Map_arPhi, 1, :)

    Map_fr = reshape(Map_fr, 1, :)
    mms = reshape(mms, 1, :)
    modes = reshape(modes, 1, :)
    Map_fiAmp = reshape(Map_fiAmp, 1, :)
    Map_fiPhi = reshape(Map_fiPhi, 1, :)
    fgridScaled = @. ifelse(
        fgrid < Map_fiAmp,
        fgrid * Map_ai + Map_bi,
        ifelse(
            fgrid < Map_fr,
            fgrid * Map_amAmp + Map_bmAmp,
            fgrid * Map_arAmp + Map_brAmp,
        ),
    )

    # Map the ampliude's range
    # We divide by the leading order l=m=2 behavior, and then scale in the expected PN behavior for the multipole of interest.

    beta_term1 = _onePointFiveSpinPN_hphc(fgrid, chi_s, chi_a, modes, mms, eta, Seta)
    beta_term2 =
        _onePointFiveSpinPN_hphc(2.0 * fgrid ./ mms, chi_s, chi_a, modes, mms, eta, Seta)
    HMamp_term1 = _onePointFiveSpinPN_hphc(fgridScaled, chi_s, chi_a, modes, mms, eta, Seta)


    HMamp_term2 = pi * sqrt(eta * 2.0 / 3.0) .* ((pi .* fgridScaled) .^ (-7.0 / 6.0))


    # The (3,3) and (4,3) modes vanish if eta=0.25 (equal mass case) and the (2,1) mode vanishes if both eta=0.25 and chi1z=chi2z
    # This results in NaNs having 0/0, correct for this using replace!()
    completeAmpl = @. Overallamp *
       amp0 *
       (fgridScaled^(-7.0 / 6.0)) *
       ifelse(
           fgridScaled < fInsJoin_Ampl,
           1.0 +
           (fgridScaled^(2.0 / 3.0)) * Acoeffs.two_thirds +
           (fgridScaled^(4.0 / 3.0)) * Acoeffs.four_thirds +
           (fgridScaled^(5.0 / 3.0)) * Acoeffs.five_thirds +
           (fgridScaled^(7.0 / 3.0)) * Acoeffs.seven_thirds +
           (fgridScaled^(8.0 / 3.0)) * Acoeffs.eight_thirds +
           fgridScaled * (
               Acoeffs.one +
               fgridScaled * Acoeffs.two +
               fgridScaled * fgridScaled * Acoeffs.three
           ),
           ifelse(
               fgridScaled < fpeak,
               delta0 +
               fgridScaled * delta1 +
               fgridScaled *
               fgridScaled *
               (delta2 + fgridScaled * delta3 + fgridScaled * fgridScaled * delta4),
               ifelse(
                   fgridScaled < fcutPar,
                   exp(-(fgridScaled - fring) * gamma2 / (fdamp * gamma3)) *
                   (fdamp * gamma3 * gamma1) / (
                       (fgridScaled - fring) * (fgridScaled - fring) +
                       fdamp * gamma3 * fdamp * gamma3
                   ),
                   0.0,
               ),
           ),
       )

    AmplsAllModes = replace!(
        completeAmpl .* (beta_term1 ./ beta_term2) .* HMamp_term1 ./ HMamp_term2,
        NaN => 0.0,
    )
    #AmplsAllModes = AmplsAllModes.transpose(0,2,1)
    #C1MRDHM, C2MRDHM, Rholm, Taulm = C1MRDHM', C2MRDHM', Rholm', Taulm'
    C1MRDHM, C2MRDHM, Rholm, Taulm = reshape(C1MRDHM, 1, :),
    reshape(C2MRDHM, 1, :),
    reshape(Rholm, 1, :),
    reshape(Taulm, 1, :)

    tmpMf = @. Map_amPhi * Map_fiPhi + Map_bmPhi
    tmpMf = reshape(tmpMf, 1, :)
    PhDBconst =
        1.0 ./ Map_amPhi .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )

    tmpMf = @. Map_arPhi * Map_fr + Map_brPhi
    tmpMf = reshape(tmpMf, 1, :)
    PhDCconst =
        1.0 ./ Map_arPhi .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )

    tmpMf = @. Map_ai * Map_fiPhi + Map_bi
    tmpMf = reshape(tmpMf, 1, :)
    PhDBAterm =
        1.0 ./ Map_ai .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )

    tmpMf = @. Map_amPhi * Map_fr + Map_bmPhi
    tmpMf = reshape(tmpMf, 1, :)
    tmpphaseC =
        -PhDBconst .+ PhDBAterm .+
        1.0 ./ Map_amPhi .* _completePhase(
            tmpMf,
            C1MRDHM,
            C2MRDHM,
            Rholm,
            Taulm,
            C1Int,
            C2Int,
            PhiInspcoeffs,
            fInsJoin_PHI,
            alpha1,
            alpha2,
            alpha3,
            alpha4,
            alpha5,
            beta1,
            beta2,
            beta3,
            etaInv,
            fMRDJoinPh,
            fcutPar,
            fring,
            fdamp,
        )



    PhisAllModes =
        ifelse.(
            fgrid .< Map_fiPhi,
            _completePhase(
                fgrid .* Map_ai .+ Map_bi,
                C1MRDHM,
                C2MRDHM,
                Rholm,
                Taulm,
                C1Int,
                C2Int,
                PhiInspcoeffs,
                fInsJoin_PHI,
                alpha1,
                alpha2,
                alpha3,
                alpha4,
                alpha5,
                beta1,
                beta2,
                beta3,
                etaInv,
                fMRDJoinPh,
                fcutPar,
                fring,
                fdamp,
            ) ./ Map_ai,
            ifelse.(
                fgrid .< Map_fr,
                -PhDBconst .+ PhDBAterm .+
                _completePhase(
                    fgrid .* Map_amPhi .+ Map_bmPhi,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin_PHI,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoinPh,
                    fcutPar,
                    fring,
                    fdamp,
                ) ./ Map_amPhi,
                -PhDCconst .+ tmpphaseC .+
                _completePhase(
                    fgrid .* Map_arPhi .+ Map_brPhi,
                    C1MRDHM,
                    C2MRDHM,
                    Rholm,
                    Taulm,
                    C1Int,
                    C2Int,
                    PhiInspcoeffs,
                    fInsJoin_PHI,
                    alpha1,
                    alpha2,
                    alpha3,
                    alpha4,
                    alpha5,
                    beta1,
                    beta2,
                    beta3,
                    etaInv,
                    fMRDJoinPh,
                    fcutPar,
                    fring,
                    fdamp,
                ) ./ Map_arPhi,
            ),
        )
    PhisAllModes = @. PhisAllModes - t0 * (fgrid - fRef) - mms * phi0 + complShiftm[mms+1]
    #modes = expand_dims(modes, len(modes.shape))
    Y, Ymstar = _spinWeighted_SphericalHarmonic(iota, modes)
    Ymstar = conj(Ymstar)
    Y, Ymstar = reshape(Y, 1, :), reshape(Ymstar, 1, :)
    ells = reshape(ells, 1, :)


    hp = vec(sum(
        AmplsAllModes .* exp.(-1im .* PhisAllModes) .*
        (0.5 .* (Y .+ ((-1) .^ ells) .* Ymstar)),
        dims = 2,
    ))  # vec() to obtain an array from a nx1 matrix
    hc = -vec(sum(
        AmplsAllModes .* exp.(-1im .* PhisAllModes) .*
        (-1im .* 0.5 .* (Y .- ((-1) .^ ells) .* Ymstar)),
        dims = 2,
    ))

    return hp, hc

end

"""
useful functions for hphc()
"""

function _onePointFiveSpinPN_Ampl(infreqs, l, m, chi_s, chi_a, eta, Seta) #is the output a number or an array?
    # PN amplitudes function, needed to scale

    v = @. (2.0 * pi * infreqs / m)^(1.0 / 3.0)
    v2 = v .* v
    v3 = v2 .* v

    if (2 == l) && (2 == m)
        Hlm = 1.0    # is it an array or a scalar? infreqs is an array, so is v

    elseif (2 == l) && (1 == m)
        Hlm = @. (sqrt(2.0) / 3.0) * (
            v * Seta - v2 * 1.5 * (chi_a + Seta * chi_s) +
            v3 * Seta * ((335.0 / 672.0) + (eta * 117.0 / 56.0)) +
            v3 *
            v *
            (
                chi_a * (3427.0 / 1344.0 - eta * 2101.0 / 336.0) +
                Seta * chi_s * (3427.0 / 1344 - eta * 965 / 336) +
                Seta * (-1im * 0.5 - pi - 2 * 1im * 0.69314718056)
            )
        )

    elseif (3 == l) && (3 == m)
        Hlm = @. 0.75 * sqrt(5.0 / 7.0) * (v * Seta)

    elseif (3 == l) && (2 == m)
        Hlm = @. (1.0 / 3.0) * sqrt(5.0 / 7.0) * (v2 * (1.0 - 3.0 * eta))

    elseif (4 == l) && (4 == m)
        Hlm = @. (4.0 / 9.0) * sqrt(10.0 / 7.0) * v2 * (1.0 - 3.0 * eta)

    elseif (4 == l) && (3 == m)
        Hlm = @. 0.75 * sqrt(3.0 / 35.0) * v3 * Seta * (1.0 - 2.0 * eta)

    else
        error("Mode not present in IMRPhenomHM waveform model.")
    end
    # Compute the final PN Amplitude at Leading Order in Mf

    return @. pi * sqrt(eta * 2.0 / 3.0) .* (v .^ (-3.5)) * abs.(Hlm)
end

"""
useful functions for hphc()
"""

function _completePhase(
    infreqs,
    C1MRDuse,
    C2MRDuse,
    RhoUse,
    TauUse,
    C1Int,
    C2Int,
    PhiInspcoeffs,
    fInsJoin_PHI,
    alpha1,
    alpha2,
    alpha3,
    alpha4,
    alpha5,
    beta1,
    beta2,
    beta3,
    etaInv,
    fMRDJoinPh,
    fcutPar,
    fring,
    fdamp,
)
    return @. ifelse(
        infreqs < fInsJoin_PHI,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (infreqs^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (infreqs^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (infreqs^(1.0 / 3.0)) * log(pi * infreqs) / 3.0 +
        PhiInspcoeffs.log * log(pi * infreqs) / 3.0 +
        PhiInspcoeffs.min_third * (infreqs^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (infreqs^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / infreqs +
        PhiInspcoeffs.min_four_thirds * (infreqs^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (infreqs^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * infreqs +
            PhiInspcoeffs.four_thirds * (infreqs^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (infreqs^(5.0 / 3.0)) +
            PhiInspcoeffs.two * infreqs * infreqs
        ) * etaInv,
        ifelse(
            infreqs < fMRDJoinPh,
            (beta1 * infreqs - beta3 / (3.0 * infreqs^3) + beta2 * log(infreqs)) * etaInv +
            C1Int +
            C2Int * infreqs,
            ifelse(
                infreqs < fcutPar,
                (
                    -(alpha2 / infreqs) +
                    (4.0 / 3.0) * (alpha3 * (infreqs^(3.0 / 4.0))) +
                    alpha1 * infreqs +
                    alpha4 *
                    RhoUse *
                    atan((infreqs - alpha5 * fring) / (fdamp * RhoUse * TauUse))
                ) * etaInv +
                C1MRDuse +
                C2MRDuse * infreqs,
                0.0,
            ),
        ),
    )
end

"""
useful functions for hphc()
"""

function _onePointFiveSpinPN_hphc(infreqs, chi_s, chi_a, modes, mms, eta, Seta) #this one is different from the one in Ampl since it computes ones for all the multipoles, to distinguish look at the number of inputs
    # PN amplitudes function, needed to scale

    v = @. (2.0 * pi * infreqs / mms)^(1.0 / 3.0)
    v2 = v .* v
    v3 = v2 .* v
    Hlm = @. ifelse(
        modes == 21,
        (sqrt(2.0) / 3.0) * (
            v * Seta - v2 * 1.5 * (chi_a + Seta * chi_s) +
            v3 * Seta * ((335.0 / 672.0) + (eta * 117.0 / 56.0)) +
            v3 *
            v *
            (
                chi_a * (3427.0 / 1344.0 - eta * 2101.0 / 336.0) +
                Seta * chi_s * (3427.0 / 1344 - eta * 965 / 336) +
                Seta * (-1im * 0.5 - pi - 2 * 1im * 0.69314718056)
            )
        ),
        ifelse(
            modes == 22,
            1.0,
            ifelse(
                modes == 32,
                (1.0 / 3.0) * sqrt(5.0 / 7.0) * (v2 * (1.0 - 3.0 * eta)),
                ifelse(
                    modes == 33,
                    0.75 * sqrt(5.0 / 7.0) * (v * Seta),
                    ifelse(
                        modes == 43,
                        0.75 * sqrt(3.0 / 35.0) * v3 * Seta * (1.0 - 2.0 * eta),
                        (4.0 / 9.0) * sqrt(10.0 / 7.0) * v2 * (1.0 - 3.0 * eta),
                    ),
                ),
            ),
        ),
    )

    # Compute the final PN Amplitude at Leading Order in Mf

    return @. pi * sqrt(eta * 2.0 / 3.0) * (v^(-3.5)) * abs(Hlm)
end

"""
useful functions for hphc()
"""

function _spinWeighted_SphericalHarmonic(theta, modes)
    # Taken from arXiv:0709.0093v3 eq. (II.7), (II.8) and LALSimulation for the s=-2 case and up to l=4.
    # We assume already phi=0 and s=-2 to simplify the function

    Ylm = @. ifelse(
        modes == 21,
        sqrt(5.0 / (16.0 * pi)) * sin(theta) * (1.0 + cos(theta)),
        ifelse(
            modes == 22,
            sqrt(5.0 / (64.0 * pi)) * (1.0 + cos(theta)) * (1.0 + cos(theta)),
            ifelse(
                modes == 32,
                sqrt(7.0 / pi) * (cos(theta * 0.5)^4) * (-2.0 + 3.0 * cos(theta)) * 0.5,
                ifelse(
                    modes == 33,
                    -sqrt(21.0 / (2.0 * pi)) * (cos(theta * 0.5)^5) * sin(theta * 0.5),
                    ifelse(
                        modes == 43,
                        -3.0 *
                        sqrt(7.0 / (2.0 * pi)) *
                        (cos(theta * 0.5)^5) *
                        (-1.0 + 2.0 * cos(theta)) *
                        sin(theta * 0.5),
                        3.0 * sqrt(7.0 / pi) * (cos(theta * 0.5)^6) * (sin(theta * 0.5)^2),
                    ),
                ),
            ),
        ),
    )
    Ylminm = @. ifelse(
        modes == 21,
        sqrt(5.0 / (16.0 * pi)) * sin(theta) * (1.0 - cos(theta)),
        ifelse(
            modes == 22,
            sqrt(5.0 / (64.0 * pi)) * (1.0 - cos(theta)) * (1.0 - cos(theta)),
            ifelse(
                modes == 32,
                sqrt(7.0 / (4.0 * pi)) *
                (2.0 + 3.0 * cos(theta)) *
                ((sin(theta * 0.5))^(4.0)),
                ifelse(
                    modes == 33,
                    sqrt(21.0 / (2.0 * pi)) * cos(theta * 0.5) * ((sin(theta * 0.5))^(5.0)),
                    ifelse(
                        modes == 43,
                        3.0 *
                        sqrt(7.0 / (2.0 * pi)) *
                        cos(theta * 0.5) *
                        (1.0 + 2.0 * cos(theta)) *
                        ((sin(theta * 0.5))^5.0),
                        3.0 *
                        sqrt(7.0 / pi) *
                        (cos(theta * 0.5) * cos(theta * 0.5)) *
                        ((sin(theta * 0.5))^6.0),
                    ),
                ),
            ),
        ),
    )

    return Ylm, Ylminm
end

"""
useful functions for hphc()
"""

function _RDfreqCalc(model::PhenomHM, finalmass, finalspin, l, m)
    """
    Compute the real and imaginary parts of the complex ringdown frequency for the :math:`(l,m)` mode as in :py:class:`LALSimIMRPhenomHM.c` line 189. This function includes all fits of the different modes.

     array or float finalmass: Mass(es) of the final object(s).
     array or float finalspin: Spin(s) of the final object(s).
     int l: :math:`l` of the chosen mode.
     int m: :math:`m` of the chosen mode.
     Real and imaginary parts of the complex ringdown frequency (ringdown and damping frequencies).
    :rtype: tuple(array, array) or tuple(float, float)

    """

    # Domain mapping for dimnesionless BH spin
    alpha = log(2.0 - finalspin) / log(3.0)
    beta = 1.0 / (2.0 + l - abs(m))
    kappa = alpha^beta
    kappa2 = kappa * kappa
    kappa3 = kappa * kappa2
    kappa4 = kappa * kappa3

    if (2 == l) && (2 == m)
        res =
            1.0 +
            kappa * (
                1.557847 * exp(2.903124 * 1im) +
                1.95097051 * exp(5.920970 * 1im) * kappa +
                2.09971716 * exp(2.760585 * 1im) * kappa2 +
                1.41094660 * exp(5.914340 * 1im) * kappa3 +
                0.41063923 * exp(2.795235 * 1im) * kappa4
            )

    elseif (3 == l) && (2 == m)
        res =
            1.022464 * exp(0.004870 * 1im) +
            0.24731213 * exp(0.665292 * 1im) * kappa +
            1.70468239 * exp(3.138283 * 1im) * kappa2 +
            0.94604882 * exp(0.163247 * 1im) * kappa3 +
            1.53189884 * exp(5.703573 * 1im) * kappa4 +
            2.28052668 * exp(2.685231 * 1im) * kappa4 * kappa +
            0.92150314 * exp(5.841704 * 1im) * kappa4 * kappa2

    elseif (4 == l) && (4 == m)
        res =
            2.0 +
            kappa * (
                2.658908 * exp(3.002787 * 1im) +
                2.97825567 * exp(6.050955 * 1im) * kappa +
                3.21842350 * exp(2.877514 * 1im) * kappa2 +
                2.12764967 * exp(5.989669 * 1im) * kappa3 +
                0.60338186 * exp(2.830031 * 1im) * kappa4
            )

    elseif (2 == l) && (1 == m)
        res =
            0.589113 * exp(0.043525 * 1im) +
            0.18896353 * exp(2.289868 * 1im) * kappa +
            1.15012965 * exp(5.810057 * 1im) * kappa2 +
            6.04585476 * exp(2.741967 * 1im) * kappa3 +
            11.12627777 * exp(5.844130 * 1im) * kappa4 +
            9.34711461 * exp(2.669372 * 1im) * kappa4 * kappa +
            3.03838318 * exp(5.791518 * 1im) * kappa4 * kappa2

    elseif (3 == l) && (3 == m)
        res =
            1.5 +
            kappa * (
                2.095657 * exp(2.964973 * 1im) +
                2.46964352 * exp(5.996734 * 1im) * kappa +
                2.66552551 * exp(2.817591 * 1im) * kappa2 +
                1.75836443 * exp(5.932693 * 1im) * kappa3 +
                0.49905688 * exp(2.781658 * 1im) * kappa4
            )

    elseif (4 == l) && (3 == m)
        res =
            1.5 +
            kappa * (
                0.205046 * exp(0.595328 * 1im) +
                3.10333396 * exp(3.016200 * 1im) * kappa +
                4.23612166 * exp(6.038842 * 1im) * kappa2 +
                3.02890198 * exp(2.826239 * 1im) * kappa3 +
                0.90843949 * exp(5.915164 * 1im) * kappa4
            )

    else
        error("Mode not present in IMRPhenomHM waveform model.")
    end

    if m < 0
        res = -conj(res)
    end
    fring = real(res) / (2.0 * pi * finalmass)

    fdamp = imag(res) / (2.0 * pi * finalmass)

    return fring, fdamp
end



##############################################################################
# IMRPhenomNSBH WAVEFORM
##############################################################################
    
"""
IMRPhenomNSBH waveform model

The inputs labelled as 1 refer to the BH (e.g. ``'chi1z'``) and with 2 to the NS (e.g. ``'Lambda2'``)

Relevant references:
    [1] `arXiv:1508.07250 <https://arxiv.org/abs/1508.07250>`_
    
    [2] `arXiv:1508.07253 <https://arxiv.org/abs/1508.07253>`_
    
    [3] `arXiv:1509.00512 <https://arxiv.org/abs/1509.00512>`_
    
    [4] `arXiv:1905.06011 <https://arxiv.org/abs/1905.06011>`_

 bool, optional verbose: Boolean specifying if the code has to print additional details during execution.
 kwargs: Optional arguments to be passed to the parent class :py:class:`WaveFormModel`, such as ``is_chi1chi2``.
    
NOTE: In LAL, to compute the parameter xi_tide in arXiv:1509.00512 eq. (8), the roots are extracted.
        In Python this would break the possibility to vectorise so, to circumvent the issue, we compute
        a grid of xi_tide as a function of the compactness, mass ratio and BH spin, and then use a 3D
        interpolator. The first time the code runs, if this interpolator is not already present, it will be
        computed (the base resolution of the grid is 200 pts per parameter, that we find
        sufficient to reproduce LAL waveforms with good precision, given the smooth behaviour of the function,
        but this can be raised if needed. In this case, it is necessary to change the name of the file assigned to self.path_xiTide_tab and the res input passed to _make_xiTide_interpolator())
"""
# All is taken from LALSimulation and arXiv:1508.07250, arXiv:1508.07253, arXiv:1509.00512, arXiv:1905.06011
    # Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase
    # Dimensionless frequency (Mf) at which we define the end of the waveform
    
""" helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Phi(model::PhenomNSBH,
    f,
    mc,
    eta,
    chi1,
    chi2, # tuned only for chi2 = 0
    Lambda,
    Lambda2;
    fInsJoin = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
)
    # if Lambda2 !==0.
    #     println("you give two values of Lambda, only one is needed for NSBH")
    #     println("the code automaticaly set the second one to zero")
    # end
    return Phi(model, f, mc, eta, chi1, chi2, Lambda; fInsJoin = fInsJoin, fcutPar = fcutPar, GMsun_over_c3 = GMsun_over_c3)
end

"""
Compute the phase of the GW as a function of frequency, given the events parameters.

    Phi(PhenomNSBH(), f, mc, eta, chi1, chi2, Lambda)

Note that the second spin is assumed to be zero.

#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the phase will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `Lambda`: Dimensionless tidal deformability of the NS.
#### Optional arguments:
-  `fInsJoin_PHI`: Dimensionless frequency (Mf) at which the inspiral phase switches to the intermediate phase. Default is 0.018.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW phase for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The phase is given in radians.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomD(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Phi(PhenomNSBH(), f, mc, eta, chi1, chi2)
```
"""


function Phi(model::PhenomNSBH,
        f,
        mc,
        eta,
        chi1,
        chi2, # tuned only for chi2 = 0
        Lambda;
        fInsJoin = 0.018,
        fcutPar = 0.2,
        GMsun_over_c3 = uc.GMsun_over_c3,
    )

    path=pwd()*"/useful_files/WFfiles/"
    QNMgrid_a = _readQNMgrid_a(path)
    QNMgrid_fring = _readQNMgrid_fring(path)
    QNMgrid_fdamp = _readQNMgrid_fdamp(path)


    M = mc / (eta^(0.6))
    eta2 = eta * eta # These can speed up a bit, we call them multiple times
    eta3 = eta * eta2
    etaInv = 1. / eta

    pi2 = pi * pi
    
    chi12, chi22 = chi1*chi1, chi2*chi2
    chi1dotchi2 = chi1*chi2

    # A non-zero tidal deformability induces a quadrupole moment (for BBH it is 1).
    # Taken from arXiv:1303.1528 eq. (54) and Tab. I
    QuadMon1 =1.
    QuadMon2 = ifelse(Lambda<1e-5, 1., exp(0.194 + 0.0936*log(Lambda) + 0.0474*log(Lambda)^2 - 0.00421*log(Lambda)^3 + 0.000123*log(Lambda)^4))
    
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)        
    SetaPlus1 = 1.0 + Seta
    chi_s = 0.5 * (chi1 + chi2)
    chi_a = 0.5 * (chi1 - chi2)
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    chi_s2, chi_a2 = chi_s*chi_s, chi_a*chi_a
    chi1dotchi2 = chi1*chi2
    chi_sdotchi_a = chi_s*chi_a
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # We work in dimensionless frequency M*f, not f
    fgrid = M*GMsun_over_c3.*f
    # As in arXiv:1508.07253 eq. (4) and LALSimIMRPhenomD_internals.c line 97
    chiPN = (chi_s * (1.0 - eta * 76.0 / 113.0) + Seta * chi_a)
    xi = - 1.0 + chiPN
    xi2 = xi * xi
    # Compute final spin and radiated energy for IMRPhenomNSBH, the rest is equivalent to IMRPhenomD_NRTidalv2
    # Get remnant spin for assumed aligned spin system, from arXiv:1903.11622 Table I and eq. (4), (5) and (6)
    
    p1_remSp = ((-5.44187381e-03*chi1 + 7.91165608e-03) + (2.33362046e-02*chi1 + 2.47764497e-02)*eta)*eta
    p2_remSp = ((-8.56844797e-07*chi1 - 2.81727682e-06) + (6.61290966e-06*chi1 + 4.28979016e-05)*eta)*eta
    p3_remSp = ((-3.04174272e-02*chi1 + 2.54889050e-01) + (1.47549350e-01*chi1 - 4.27905832e-01)*eta)*eta
    
    modelRemSp = (1. + Lambda * p1_remSp + Lambda^2 * p2_remSp) / ((1. + Lambda*p3_remSp^2)*(1. + Lambda*p3_remSp^2))

    modelRemSp = ifelse((chi1 < 0.) & (eta < 0.188), 1., modelRemSp)
    modelRemSp = ifelse(chi1 < -0.5, 1., modelRemSp)
    modelRemSp = ifelse(modelRemSp > 1., 1., modelRemSp)
    
    #del p1_remSp, p2_remSp, p3_remSp
    
    # Work with spin variables weighted on square of the BH mass over total mass
    S1BH = chi1 * m1ByM^2
    Shat = S1BH / (m1ByM^2 + m2ByM^2) # this would be = (chi1*m1*m1 + chi2*m2*m2)/(m1*m1 + m2*m2), but chi2=0 by assumption
    
    # Compute fit to L_orb in arXiv:1611.00332 eq. (16)
    Lorb = (2. *sqrt(3.)*eta + 5.24*3.8326341618708577*eta2 + 1.3*(-9.487364155598392)*eta3)/(1. + 2.88*2.5134875145648374*eta) + ((-0.194)*1.0009563702914628*Shat*(4.409160174224525*eta + 0.5118334706832706*eta2 + (64. - 16. *4.409160174224525 - 4. *0.5118334706832706)*eta3) 
    + 0.0851*0.7877509372255369*Shat^2*(8.77367320110712*eta + (-32.060648277652994)*eta2 + (64. - 16. *8.77367320110712 - 4. *(-32.060648277652994))*eta3) + 0.00954*0.6540138407185817*Shat^3*(22.830033250479833*eta + (-153.83722669033995)*eta2 + (64. - 16. *22.830033250479833 - 4. *(-153.83722669033995))*eta3))/(1. + (-0.579)*0.8396665722805308*Shat*(1.8804718791591157 + (-4.770246856212403)*eta + 0. *eta2 + (64. - 64. *1.8804718791591157 - 16. *(-4.770246856212403) - 4. *0.)*eta3)) + 0.3223660562764661*Seta*eta2*(1. + 9.332575956437443*eta)*chi1 + 2.3170397514509933*Shat*Seta*eta3*(1. + (-3.2624649875884852)*eta)*chi1 + (-0.059808322561702126)*eta3*chi12;
    
    chif = (Lorb + S1BH)*modelRemSp

    Erad = _radiatednrg(model, eta, chi1, chi2)
    # Compute ringdown and damping frequencies from interpolators
    fring = LinearInterpolation(QNMgrid_a, QNMgrid_fring)(real(chif)) / (1.0 - Erad)
    fdamp = LinearInterpolation(QNMgrid_a, QNMgrid_fdamp)(real(chif)) / (1.0 - Erad)

    # Compute sigma coefficients appearing in arXiv:1508.07253 eq. (28)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    sigma1 =
        2096.551999295543 +
        1463.7493168261553 * eta +
        (
            1312.5493286098522 + 18307.330017082117 * eta - 43534.1440746107 * eta2 +
            (-833.2889543511114 + 32047.31997183187 * eta - 108609.45037520859 * eta2) *
            xi +
            (452.25136398112204 + 8353.439546391714 * eta - 44531.3250037322 * eta2) * xi2
        ) * xi
    sigma2 =
        -10114.056472621156 - 44631.01109458185 * eta +
        (
            -6541.308761668722 - 266959.23419307504 * eta +
            686328.3229317984 * eta2 +
            (3405.6372187679685 - 437507.7208209015 * eta + 1.6318171307344697e6 * eta2) *
            xi +
            (-7462.648563007646 - 114585.25177153319 * eta + 674402.4689098676 * eta2) * xi2
        ) * xi
    sigma3 =
        22933.658273436497 +
        230960.00814979506 * eta +
        (
            14961.083974183695 + 1.1940181342318142e6 * eta - 3.1042239693052764e6 * eta2 +
            (-3038.166617199259 + 1.8720322849093592e6 * eta - 7.309145012085539e6 * eta2) *
            xi +
            (42738.22871475411 + 467502.018616601 * eta - 3.064853498512499e6 * eta2) * xi2
        ) * xi
    sigma4 =
        -14621.71522218357 - 377812.8579387104 * eta +
        (
            -9608.682631509726 - 1.7108925257214056e6 * eta +
            4.332924601416521e6 * eta2 +
            (
                -22366.683262266528 - 2.5019716386377467e6 * eta +
                1.0274495902259542e7 * eta2
            ) * xi +
            (-85360.30079034246 - 570025.3441737515 * eta + 4.396844346849777e6 * eta2) *
            xi2
        ) * xi

    # Compute beta coefficients appearing in arXiv:1508.07253 eq. (16)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    beta1 =
        97.89747327985583 - 42.659730877489224 * eta +
        (
            153.48421037904913 - 1417.0620760768954 * eta +
            2752.8614143665027 * eta2 +
            (138.7406469558649 - 1433.6585075135881 * eta + 2857.7418952430758 * eta2) *
            xi +
            (41.025109467376126 - 423.680737974639 * eta + 850.3594335657173 * eta2) * xi2
        ) * xi
    beta2 =
        -3.282701958759534 - 9.051384468245866 * eta +
        (
            -12.415449742258042 + 55.4716447709787 * eta - 106.05109938966335 * eta2 +
            (-11.953044553690658 + 76.80704618365418 * eta - 155.33172948098394 * eta2) *
            xi +
            (-3.4129261592393263 + 25.572377569952536 * eta - 54.408036707740465 * eta2) *
            xi2
        ) * xi
    beta3 =
        -0.000025156429818799565 +
        0.000019750256942201327 * eta +
        (
            -0.000018370671469295915 +
            0.000021886317041311973 * eta +
            0.00008250240316860033 * eta2 +
            (
                7.157371250566708e-6 - 0.000055780000112270685 * eta +
                0.00019142082884072178 * eta2
            ) * xi +
            (
                5.447166261464217e-6 - 0.00003220610095021982 * eta +
                0.00007974016714984341 * eta2
            ) * xi2
        ) * xi

    # Compute alpha coefficients appearing in arXiv:1508.07253 eq. (14)
    # They derive from a fit, whose numerical coefficients are in arXiv:1508.07253 Tab. 5
    alpha1 =
        43.31514709695348 +
        638.6332679188081 * eta +
        (
            -32.85768747216059 + 2415.8938269370315 * eta - 5766.875169379177 * eta2 +
            (-61.85459307173841 + 2953.967762459948 * eta - 8986.29057591497 * eta2) * xi +
            (-21.571435779762044 + 981.2158224673428 * eta - 3239.5664895930286 * eta2) *
            xi2
        ) * xi
    alpha2 =
        -0.07020209449091723 - 0.16269798450687084 * eta +
        (
            -0.1872514685185499 + 1.138313650449945 * eta - 2.8334196304430046 * eta2 +
            (-0.17137955686840617 + 1.7197549338119527 * eta - 4.539717148261272 * eta2) *
            xi +
            (-0.049983437357548705 + 0.6062072055948309 * eta - 1.682769616644546 * eta2) *
            xi2
        ) * xi
    alpha3 =
        9.5988072383479 - 397.05438595557433 * eta +
        (
            16.202126189517813 - 1574.8286986717037 * eta +
            3600.3410843831093 * eta2 +
            (27.092429659075467 - 1786.482357315139 * eta + 5152.919378666511 * eta2) * xi +
            (11.175710130033895 - 577.7999423177481 * eta + 1808.730762932043 * eta2) * xi2
        ) * xi
    alpha4 =
        -0.02989487384493607 +
        1.4022106448583738 * eta +
        (
            -0.07356049468633846 +
            0.8337006542278661 * eta +
            0.2240008282397391 * eta2 +
            (-0.055202870001177226 + 0.5667186343606578 * eta + 0.7186931973380503 * eta2) *
            xi +
            (
                -0.015507437354325743 +
                0.15750322779277187 * eta +
                0.21076815715176228 * eta2
            ) * xi2
        ) * xi
    alpha5 =
        0.9974408278363099 - 0.007884449714907203 * eta +
        (
            -0.059046901195591035 + 1.3958712396764088 * eta - 4.516631601676276 * eta2 +
            (-0.05585343136869692 + 1.7516580039343603 * eta - 5.990208965347804 * eta2) *
            xi +
            (-0.017945336522161195 + 0.5965097794825992 * eta - 2.0608879367971804 * eta2) *
            xi2
        ) * xi

    # Compute the TF2 phase coefficients and put them in a dictionary (spin effects are included up to 3.5PN)
    TF2OverallAmpl = 3 / (128.0 * eta)

    # For 3PN coeff we use chi1 and chi2 so to have the quadrupole moment explicitly appearing

    TF2_5coeff_tmp =
        38645.0 * pi / 756.0 - 65.0 * pi * eta / 9.0 - (
            (732985.0 / 2268.0 - 24260.0 * eta / 81.0 - 340.0 * eta2 / 9.0) * chi_s +
            (732985.0 / 2268.0 + 140.0 * eta / 9.0) * Seta * chi_a
        ) #variable to be used later
    TF2_6coeff_tmp =
        11583.231236531 / 4.694215680 - 640.0 / 3.0 * pi2 -
        684.8 / 2.1 * MathConstants.eulergamma +
        eta * (-15737.765635 / 3.048192 + 225.5 / 1.2 * pi2) +
        eta2 * 76.055 / 1.728 - eta2 * eta * 127.825 / 1.296 - log(4.0) * 684.8 / 2.1 +
        pi * chi1 * m1ByM * (1490.0 / 3.0 + m1ByM * 260.0) +
        pi * chi2 * m2ByM * (1490.0 / 3.0 + m2ByM * 260.0) +
        (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) *
        m1ByM^2*
        QuadMon1 *
        chi12 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 ) *
        m1ByM^2 *
        chi12 +
        (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) *
        m2ByM^2 *
        QuadMon2 *
        chi22 +
        (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 ) *
        m2ByM^2 *
        chi22

    TF2coeffs = TF2coeffsStructure(
        1.0,
        0.0,
        3715.0 / 756.0 + (55.0 * eta) / 9.0,
        -16.0 * pi +
        (113.0 * Seta * chi_a) / 3.0 +
        (113.0 / 3.0 - (76.0 * eta) / 3.0) * chi_s,
        5.0 * (3058.673 / 7.056 + 5429.0 / 7.0 * eta + 617.0 * eta2) / 72.0 +
        247.0 / 4.8 * eta * chi1dotchi2 - 721.0 / 4.8 * eta * chi1dotchi2 +
        (-720.0 / 9.6 * QuadMon1 + 1.0 / 9.6) * m1ByM^2  * chi12 +
        (-720.0 / 9.6 * QuadMon2 + 1.0 / 9.6) * m2ByM^2  * chi22 +
        (240.0 / 9.6 * QuadMon1 - 7.0 / 9.6) * m1ByM^2  * chi12 +
        (240.0 / 9.6 * QuadMon2 - 7.0 / 9.6) * m2ByM^2  * chi22,
        TF2_5coeff_tmp,
        TF2_5coeff_tmp * 3.0,
        TF2_6coeff_tmp - (
            (326.75 / 1.12 + 557.5 / 1.8 * eta) * eta * chi1dotchi2 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m1ByM - 120.0 * m1ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m1ByM + 125.5 / 3.6 * m1ByM^2 )
            ) *
            m1ByM *
            m1ByM *
            chi12 +
            (
                (4703.5 / 8.4 + 2935.0 / 6.0 * m2ByM - 120.0 * m2ByM^2 ) +
                (-4108.25 / 6.72 - 108.5 / 1.2 * m2ByM + 125.5 / 3.6 * m2ByM^2 )
            ) *
            m2ByM *
            m2ByM *
            chi22
        ),
        -6848.0 / 21.0,
        77096675.0 * pi / 254016.0 + 378515.0 * pi * eta / 1512.0 -
        74045.0 * pi * eta2 / 756.0 +
        (
            -25150083775.0 / 3048192.0 + 10566655595.0 * eta / 762048.0 -
            1042165.0 * eta2 / 3024.0 + 5345.0 * eta2 * eta / 36.0
        ) * chi_s +
        Seta * (
            (
                -25150083775.0 / 3048192.0 + 26804935.0 * eta / 6048.0 -
                1985.0 * eta2 / 48.0
            ) * chi_a
        ),
    )
    
    PhiInspcoeffs = PhiInspcoeffsStructure(
        TF2coeffs.five * TF2OverallAmpl,
        TF2coeffs.seven * TF2OverallAmpl * (pi^(2.0 / 3.0)),
        TF2coeffs.six * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.six_log * TF2OverallAmpl * (pi^(1.0 / 3.0)),
        TF2coeffs.five_log * TF2OverallAmpl,
        TF2coeffs.four * TF2OverallAmpl * (pi^(-1.0 / 3.0)),
        TF2coeffs.three * TF2OverallAmpl * (pi^(-2.0 / 3.0)),
        TF2coeffs.two * TF2OverallAmpl / pi,
        TF2coeffs.one * TF2OverallAmpl * (pi^(-4.0 / 3.0)),
        TF2coeffs.zero * TF2OverallAmpl * (pi^(-5.0 / 3.0)),
        sigma1,
        sigma2 * 0.75,
        sigma3 * 0.6,
        sigma4 * 0.5,
    )
    
    #Now compute the coefficients to align the three parts
    
    fMRDJoin = 0.5*fring
    
    # First the Inspiral - Intermediate: we compute C1Int and C2Int coeffs
    # Equations to solve for to get C(1) continuous join
    # PhiIns (f)  =   PhiInt (f) + C1Int + C2Int f
    # Joining at fInsJoin
    # PhiIns (fInsJoin)  =   PhiInt (fInsJoin) + C1Int + C2Int fInsJoin
    # PhiIns'(fInsJoin)  =   PhiInt'(fInsJoin) + C2Int
    # This is the first derivative wrt f of the inspiral phase computed at fInsJoin, first add the PN contribution and then the higher order calibrated terms
    DPhiIns =
        (
            2.0 * TF2coeffs.seven * TF2OverallAmpl * ((pi * fInsJoin)^(7.0 / 3.0)) +
            (
                TF2coeffs.six * TF2OverallAmpl +
                TF2coeffs.six_log * TF2OverallAmpl * (1.0 + log(pi * fInsJoin) / 3.0)
            ) * ((pi * fInsJoin)^(2.0)) +
            TF2coeffs.five_log * TF2OverallAmpl * ((pi * fInsJoin)^(5.0 / 3.0)) -
            TF2coeffs.four * TF2OverallAmpl * ((pi * fInsJoin)^(4.0 / 3.0)) -
            2.0 * TF2coeffs.three * TF2OverallAmpl * (pi * fInsJoin) -
            3.0 * TF2coeffs.two * TF2OverallAmpl * ((pi * fInsJoin)^(2.0 / 3.0)) -
            4.0 * TF2coeffs.one * TF2OverallAmpl * ((pi * fInsJoin)^(1.0 / 3.0)) -
            5.0 * TF2coeffs.zero * TF2OverallAmpl
        ) * pi / (3.0 * ((pi * fInsJoin)^(8.0 / 3.0))) +
        (
            sigma1 +
            sigma2 * (fInsJoin^(1.0 / 3.0)) +
            sigma3 * (fInsJoin^(2.0 / 3.0)) +
            sigma4 * fInsJoin
        ) * etaInv

    # This is the first derivative of the Intermediate phase computed at fInsJoin
    DPhiInt = (beta1 + beta3 / (fInsJoin^4) + beta2 / fInsJoin) * etaInv
    C2Int = DPhiIns - DPhiInt

    # This is the inspiral phase computed at fInsJoin
    PhiInsJoin =
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * (fInsJoin^(2.0 / 3.0)) +
        PhiInspcoeffs.third * (fInsJoin^(1.0 / 3.0)) +
        PhiInspcoeffs.third_log * (fInsJoin^(1.0 / 3.0)) * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.log * log(pi * fInsJoin) / 3.0 +
        PhiInspcoeffs.min_third * (fInsJoin^(-1.0 / 3.0)) +
        PhiInspcoeffs.min_two_thirds * (fInsJoin^(-2.0 / 3.0)) +
        PhiInspcoeffs.min_one / fInsJoin +
        PhiInspcoeffs.min_four_thirds * (fInsJoin^(-4.0 / 3.0)) +
        PhiInspcoeffs.min_five_thirds * (fInsJoin^(-5.0 / 3.0)) +
        (
            PhiInspcoeffs.one * fInsJoin +
            PhiInspcoeffs.four_thirds * (fInsJoin^(4.0 / 3.0)) +
            PhiInspcoeffs.five_thirds * (fInsJoin^(5.0 / 3.0)) +
            PhiInspcoeffs.two * fInsJoin * fInsJoin
        ) * etaInv
    # This is the Intermediate phase computed at fInsJoin
    PhiIntJoin =
        beta1 * fInsJoin - beta3 / (3.0 * fInsJoin * fInsJoin * fInsJoin) +
        beta2 * log(fInsJoin)

    C1Int = PhiInsJoin - PhiIntJoin * etaInv - C2Int * fInsJoin

    # Now the same for Intermediate - Merger-Ringdown: we also need a temporary Intermediate Phase function
    PhiIntTempVal =
        (beta1 * fMRDJoin - beta3 / (3.0 * fMRDJoin^3) + beta2 * log(fMRDJoin)) * etaInv +
        C1Int +
        C2Int * fMRDJoin
    DPhiIntTempVal = C2Int + (beta1 + beta3 / (fMRDJoin^4) + beta2 / fMRDJoin) * etaInv
    DPhiMRDVal =
        (
            alpha1 +
            alpha2 / (fMRDJoin^2) +
            alpha3 / (fMRDJoin^(0.25)) +
            alpha4 / (
                fdamp * (
                    1.0 +
                    (fMRDJoin - alpha5 * fring) * (fMRDJoin - alpha5 * fring) / (fdamp^2)
                )
            )
        ) * etaInv
    PhiMRJoinTemp =
        -(alpha2 / fMRDJoin) +
        (4.0 / 3.0) * (alpha3 * (fMRDJoin^(0.75))) +
        alpha1 * fMRDJoin +
        alpha4 * atan((fMRDJoin - alpha5 * fring) / fdamp)

    C2MRD = DPhiIntTempVal - DPhiMRDVal
    C1MRD = PhiIntTempVal - PhiMRJoinTemp * etaInv - C2MRD * fMRDJoin

    fpeak = fgrid[end] # In LAL the maximum of the grid is used to rescale
    
    t0 = (alpha1 + alpha2/(fpeak*fpeak) + alpha3/(fpeak^(1. /4.)) + alpha4/(fdamp*(1. + (fpeak - alpha5*fring)*(fpeak - alpha5*fring)/(fdamp*fdamp))))/eta
    
    # LAL sets fRef as the minimum frequency, do the same
    fRef   = fgrid[1]

    phiRef = ifelse(
        fRef < fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fRef^(2.0 / 3.0) +
        PhiInspcoeffs.third * fRef^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fRef^(1.0 / 3.0) * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.log * log(pi * fRef) / 3.0 +
        PhiInspcoeffs.min_third * fRef^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fRef^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fRef +
        PhiInspcoeffs.min_four_thirds * fRef^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fRef^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fRef +
            PhiInspcoeffs.four_thirds * fRef^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fRef^(5.0 / 3.0) +
            PhiInspcoeffs.two * fRef * fRef
        ) * etaInv,
        ifelse(
            fRef < fMRDJoin,
            (beta1 * fRef - beta3 / (3.0 * fRef * fRef * fRef) + beta2 * log(fRef)) *
            etaInv +
            C1Int +
            C2Int * fRef,
            ifelse(
                fRef < fcutPar,
                (
                    -(alpha2 / fRef) +
                    (4.0 / 3.0) * (alpha3 * (fRef^(3.0 / 4.0))) +
                    alpha1 * fRef +
                    alpha4 * atan((fRef - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fRef,
                0.0,
            ),
        ),
    )

    phis = @. ifelse.(
        fgrid .< fInsJoin,
        PhiInspcoeffs.initial_phasing +
        PhiInspcoeffs.two_thirds * fgrid^(2.0 / 3.0) +
        PhiInspcoeffs.third * fgrid^(1.0 / 3.0) +
        PhiInspcoeffs.third_log * fgrid^(1.0 / 3.0) * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.log * log(pi * fgrid) / 3.0 +
        PhiInspcoeffs.min_third * fgrid^(-1.0 / 3.0) +
        PhiInspcoeffs.min_two_thirds * fgrid^(-2.0 / 3.0) +
        PhiInspcoeffs.min_one / fgrid +
        PhiInspcoeffs.min_four_thirds * fgrid^(-4.0 / 3.0) +
        PhiInspcoeffs.min_five_thirds * fgrid^(-5.0 / 3.0) +
        (
            PhiInspcoeffs.one * fgrid +
            PhiInspcoeffs.four_thirds * fgrid^(4.0 / 3.0) +
            PhiInspcoeffs.five_thirds * fgrid^(5.0 / 3.0) +
            PhiInspcoeffs.two * fgrid * fgrid
        ) * etaInv,
        ifelse(
            fgrid < fMRDJoin,
            (beta1 * fgrid - beta3 / (3.0 * fgrid * fgrid * fgrid) + beta2 * log(fgrid)) * etaInv +
            C1Int +
            C2Int * fgrid,
            ifelse(
                fgrid < fcutPar,
                (
                    -(alpha2 / fgrid) +
                    (4.0 / 3.0) * (alpha3 * (fgrid^(3.0 / 4.0))) +
                    alpha1 * fgrid +
                    alpha4 * atan((fgrid - alpha5 * fring) / fdamp)
                ) * etaInv +
                C1MRD +
                C2MRD * fgrid,
                0.0,
            ),
        ),
    )
    
    # Add the tidal contribution to the phase, as in arXiv:1905.06011
    # Compute the tidal coupling constant, arXiv:1905.06011 eq. (8) using Lambda = 2/3 k_2/C^5 (eq. (10))

    kappa2T = (3.0/13.0) * ((1.0 + 12.0*m1ByM/m2ByM)*(m2ByM^5)*Lambda)
    
    c_Newt   = 2.4375
    n_1      = -12.615214237993088
    n_3over2 =  19.0537346970349
    n_2      = -21.166863146081035
    n_5over2 =  90.55082156324926
    n_3      = -60.25357801943598
    d_1      = -15.11120782773667
    d_3over2 =  22.195327350624694
    d_2      =   8.064109635305156

    numTidal = @. 1.0 + (n_1 * ((pi*fgrid)^(2. /3.))) + (n_3over2 * pi*fgrid) + (n_2 * ((pi*fgrid)^(4. /3.))) + (n_5over2 * ((pi*fgrid)^(5. /3.))) + (n_3 * (pi*fgrid)^2)
    denTidal = @. 1.0 + (d_1 * ((pi*fgrid)^(2. /3.))) + (d_3over2 * pi*fgrid) + (d_2 * ((pi*fgrid)^(4. /3.)))
    
    tidal_phase = @. - kappa2T * c_Newt / (m1ByM * m2ByM) * ((pi*fgrid)^(5. /3.)) * numTidal / denTidal
    # println("phis: ", phis)
    # println("tidal_phase ", tidal_phase)
    # println("t0*(fgrid - fRef) ", t0.*(fgrid .- fRef))
    # println("phiRef ", phiRef)
    # println("t0 ", t0)
    # println("fgrid ", fgrid)
    # println("fRef ", fRef)
    # This pi factor is needed to include LAL fRef rescaling, so to end up with the exact same waveform
    return @. phis + ifelse(fgrid < fcutPar, - t0*(fgrid - fRef) - phiRef + pi +  tidal_phase, 0.)
end

""" helper function to do function overloading (i.e., to have different functions with the same name but different input arguments) 
"""
function Ampl(model::PhenomNSBH,
    f,
    mc,
    eta,
    chi1,
    chi2, #assumed to be chi2=0
    dL,
    Lambda1,
    Lambda2;
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)
    # if Lambda2 !==0.
    #     println("you give two values of Lambda, only one is needed for NSBH")
    #     println("the code automaticaly set the second one to zero")
    # end
    return Ampl(model, f, mc, eta, chi1, chi2, dL, Lambda1, fcutPar=fcutPar, GMsun_over_c3=GMsun_over_c3, GMsun_over_c2_Gpc=GMsun_over_c2_Gpc)
end

"""
Compute the amplitude of the GW as a function of frequency, given the events parameters.

    Ampl(PhenomNSBH(), f, mc, eta,  chi1, chi2, Lambda)
Note that the second spin is assumed to be zero.
#### Input arguments:
-  `model`: Model type, it indicates the waveform model to be used.
-  `f`: Frequency on which the amplitude will be computed, in Hz. The function accepts both a single value or an array.
-  `mc`: Chirp mass of the binary, in solar masses.
-  `eta`: Symmetric mass ratio of the binary.
-  `chi1`: Dimensionless spin of the first BH.
-  `chi2`: Dimensionless spin of the second BH.
-  `Lambda`: Dimensionless tidal deformability of the NS.
#### Optional arguments:
-  `fInsJoin_Ampl`: Dimensionless frequency (Mf) at which the inspiral amplitude switches to the intermediate amplitude. Default is 0.014.
-  `fcutPar`: Dimensionless frequency (Mf) at which we define the end of the waveform. Default is 0.2. 
#### Return:
-  GW amplitude for the chosen events evaluated for that frequency. The function returns an array of the same length as the input frequency array, if `f` is a Float, return a Float. The amplitude is dimensionless.

#### Example:
```julia
    mc = 30.
    eta = 0.25
    dL = 8.
    chi1 = 0.5
    chi2 = 0.5
    fcut = _fcut(PhenomD(), mc, eta)
    f = 10 .^(range(1.,log10(fcut),length=1000))
    Ampl(PhenomNSBH(), f, mc, eta, chi1, chi2, dL, Lambda)
```
"""


function Ampl(model::PhenomNSBH,
    f,
    mc,
    eta,
    chi1,
    chi2, #assumed to be chi2=0
    dL,
    Lambda;
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
)
    """
    Compute the amplitude of the GW as a function of frequency, given the events parameters.
    
     array f: Frequency grid on which the phase will be computed, in :math:`\\rm Hz`.
     dict(array, array, ...) kwargs: Dictionary with arrays containing the parameters of the events to compute the amplitude of, as in :py:data:`events`.
     GW amplitude for the chosen events evaluated on the frequency grid.
    :rtype: array
    
    """
    # Useful quantities

    M = mc / (eta^(3.0 / 5.0))
    eta2 = eta * eta # This can speed up a bit, we call it multiple times
    eta3 = eta * eta2
    pi2 = pi * pi
    chi12, chi22 = chi1*chi1, chi2*chi2
    
    Seta = ifelse(eta<.25,sqrt(1.0 - 4.0 * eta),0.)
    q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    # We work in dimensionless frequency M*f, not f
    fgrid = M * GMsun_over_c3 .* f
    # These are m1/Mtot and m2/Mtot
    m1ByM = 0.5 * (1.0 + Seta)
    m2ByM = 0.5 * (1.0 - Seta)
    # As in arXiv:0909.2867
    chieff = m1ByM * chi1 + m2ByM * chi2
    chisum = 2. *chieff
    chiprod = chieff^2
    
    # compute needed IMRPhenomC attributes
    # First the SPA part, LALSimIMRPhenomC_internals.c line 38
    # Frequency-domain Amplitude coefficients
    xdotaN = 64. *eta/5.
    xdota2 = -7.43/3.36 - 11. *eta/4.
    xdota3 = 4. *pi - 11.3*chieff/1.2 + 19. *eta*chisum/6.
    xdota4 = 3.4103/1.8144 + 5*chiprod + eta*(13.661/2.016 - chiprod/8.) + 5.9*eta2/1.8
    xdota5 = -pi*(41.59/6.72 + 189. *eta/8.) - chieff*(31.571/1.008 - 116.5*eta/2.4) + chisum*(21.863*eta/1.008 - 79. *eta2/6.) - 3*chieff*chiprod/4. + 9. *eta*chieff*chiprod/4.
    xdota6 = 164.47322263/1.39708800 - 17.12*MathConstants.eulergamma/1.05 + 16. *pi*pi/3 - 8.56*log(16.)/1.05 + eta*(45.1*pi*pi/4.8 - 561.98689/2.17728) + 5.41*eta2/8.96 - 5.605*eta3/2.592 - 80. *pi*chieff/3. + eta*chisum*(20. *pi/3. - 113.5*chieff/3.6) + chiprod*(64.153/1.008 - 45.7*eta/3.6) - chiprod*(7.87*eta/1.44 - 30.37*eta2/1.44)
    xdota6log = -856. /105.
    xdota7 = -pi*(4.415/4.032 - 358.675*eta/6.048 - 91.495*eta2/1.512) - chieff*(252.9407/2.7216 - 845.827*eta/6.048 + 415.51*eta2/8.64) + chisum*(158.0239*eta/5.4432 - 451.597*eta2/6.048 + 20.45*eta3/4.32 + 107. *eta*chiprod/6. - 5. *eta2*chiprod/24.) + 12. *pi*chiprod - chiprod*chieff*(150.5/2.4 + eta/8.) + chieff*chiprod*(10.1*eta/2.4 + 3. *eta2/8.)
    # Time-domain amplitude coefficients, which also enters the fourier amplitude in this model
    AN = 8. *eta*sqrt(pi/5.)
    A2 = (-107. + 55. *eta)/42.
    A3 = 2. *pi - 4. *chieff/3. + 2. *eta*chisum/3.
    A4 = -2.173/1.512 - eta*(10.69/2.16 - 2. *chiprod) + 2.047*eta2/1.512
    A5 = -10.7*pi/2.1 + eta*(3.4*pi/2.1)
    A5imag = -24. *eta
    A6 = 270.27409/6.46800 - 8.56*MathConstants.eulergamma/1.05 + 2. *pi*pi/3. + eta*(4.1*pi*pi/9.6 - 27.8185/3.3264) - 20.261*eta2/2.772 + 11.4635*eta3/9.9792 - 4.28*log(16.)/1.05
    A6log = -428. /105.
    A6imag = 4.28*pi/1.05
    
    z701, z702, z711, z710, z720 = 4.149e+00, -4.070e+00, -8.752e+01, -4.897e+01, 6.665e+02
    z801, z802, z811, z810, z820 = -5.472e-02, 2.094e-02, 3.554e-01, 1.151e-01, 9.640e-01
    z901, z902, z911, z910, z920 = -1.235e+00, 3.423e-01, 6.062e+00, 5.949e+00, -1.069e+01
    
    g1 = z701 * chieff + z702 * chiprod + z711 * eta * chieff + z710 * eta + z720 * eta2
    g1 = ifelse(g1 < 0., 0., g1)
    
    del1 = z801 * chieff + z802 * chiprod + z811 * eta * chieff + z810 * eta + z820 * eta2
    del2 = z901 * chieff + z902 * chiprod + z911 * eta * chieff + z910 * eta + z920 * eta2
    del1 = ifelse(del1 < 0., 0., del1)
    del2 = ifelse(del2 < 1.0e-4, 1.0e-4, del2)
    
    d0 = 0.015
    
    # All the other coefficients from IMRPhenomC are not needed
    
    # Now compute NSBH coefficients
    # Get NS compactness and baryonic mass, see arXiv:1608.02582 eq. (78)
    a0Comp = 0.360
    a1Comp = -0.0355
    a2Comp = 0.000705
    
    comp = ifelse(Lambda > 1., a0Comp + a1Comp*log(Lambda) + a2Comp*log(Lambda)^2, 0.5 + (3. *a0Comp-a1Comp-1.5)*Lambda^2 + (-2. *a0Comp+a1Comp+1.)*Lambda^3)
    
    # Get baryonic mass of the torus remnant of a BH-NS merger in units of the NS baryonic mass,
    # see arXiv:1509.00512 eq. (11)
    alphaTor = 0.296
    betaTor = 0.171
    # In LAL the relation is inverted each time, but this would break the vectorisation,
    # we use an interpolator on a grid of Comp, q, chi instead. Already with 100 pts per parameter the
    # agreement we find with LAL waveforms is at machine precision
    
    #xiTide = _xiTide_solver(comp, q, chi1)

    # xiTide = find_zero( (x,p)-> x^10 - 3*p[2]*x^8 + 2. *p[3]*(p[2]^(3. /2.)) *x^7 - 3. *p[1]*x^4 + 6. *p[1]*p[2]*x^2 - 3. *p[1]*p[2]^2*p[3]^2, 20., Order1(),[q, q*comp, chi1], maxiters=200)^2
    # println("xiTide: ", xiTide)
    # println(size(xiTide))
    # println(typeof(xiTide))
    #println((xiTide.value).value)
    # println((xiTide[1].partials)[1].value)
    # println((xiTide[1].partials)[2].value)
    # println((xiTide[1].partials)[3].value)
    # println((xiTide[1].partials)[4].value)



    # function zeross(eta, Lambda, chi)
    #     #println(Lambda)
    #     comp = ifelse(Lambda > 1., a0Comp + a1Comp*log(Lambda) + a2Comp*log(Lambda)^2, 0.5 + (3. *a0Comp-a1Comp-1.5)*Lambda^2 + (-2. *a0Comp+a1Comp+1.)*Lambda^3)
    #     #println(comp)
    #     Seta = ifelse(eta>=.25,0.,sqrt(1.0 - 4.0 * eta))
    #     q = 0.5*(1.0 + Seta - 2.0*eta)/eta
    #     #println(q)
    #     #println(q.partials[1])
    #     mu = comp*q
    #    c= find_zero( (x,p)-> x^10 - 3*p[2]*x^8 + 2. *p[3]*(p[2]^(3. /2.)) *x^7 - 3. *p[1]*x^4 + 6. *p[1]*p[2]*x^2 - 3. *p[1]*p[2]^2*p[3]^2
    #                        ,20., Order1(),[q, mu, chi], maxiters=100)^2
    #     partials = [(c.partials)[i].value for i in eachindex(c.partials)]
    #     cc = ForwardDiff.Dual{typeof(chi).parameters[1]}((c.value).value, partials...)
    #     return cc
    # end


    # function zeross(q,mu, chi)
    #     # println(q)
    #     # println(q.partials[1])
        

    #         println(partials...)
    #         return cc
    # end

    if typeof(eta) == Float64
        xiTide = find_zero( (x,p)-> x^10 - 3*p[2]*x^8 + 2. *p[3]*(p[2]^(3. /2.)) *x^7 - 3. *p[1]*x^4 + 6. *p[1]*p[2]*x^2 - 3. *p[1]*p[2]^2*p[3]^2, 20., Order1(),[q, q*comp, chi1], maxiters=100)^2
    else    # else we are doing the Fisher matrix, i.e. we are using FowardDiff.Duals
        #partials = [(xiTide.partials)[i].value for i in eachindex(xiTide.partials)]
        #xiTide = ForwardDiff.Dual{typeof(eta).parameters[1]}((xiTide.value).value, partials...)
        #xiTide = zeross(q, q*comp, chi1)    # rewrite w/o function

        # if we are doing the Fisher matrix we need to derive the roots of the following equation w.r.t. q, q*comp, chi1
        # to do that we need to re-pack the Duals in a way that ForwardDiff can handle
        tmp = find_zero( (x,p)-> x^10 - 3*p[2]*x^8 + 2. *p[3]*(p[2]^(3. /2.)) *x^7 - 3. *p[1]*x^4 + 6. *p[1]*p[2]*x^2 - 3. *p[1]*p[2]^2*p[3]^2
                        ,20., Order1(),[q, q*comp, chi1], maxiters=100)^2
        partials = [(tmp.partials)[i].value for i in eachindex(tmp.partials)]
        xiTide = ForwardDiff.Dual{typeof(chi1).parameters[1]}((tmp.value).value, partials...)



    end
        #xiTide = ForwardDiff.Dual{typeof(eta).parameters[1]}((xiTide[1].value).value, (xiTide[1].partials)[1].value, (xiTide[1].partials)[2].value, (xiTide[1].partials)[3].value, (xiTide[1].partials)[4].value, (xiTide[1].partials)[5].value, (xiTide[1].partials)[6].value)#, (xiTide[1].partials)[7].value, (xiTide[1].partials)[8].value, (xiTide[1].partials)[9].value, (xiTide[1].partials)[10].value, (xiTide[1].partials)[11].value)
    #, xiTide[2].value, xiTide[3].value, xiTide[4].value, xiTide[5].value, xiTide[6].value, xiTide[7].value, xiTide[8].value, xiTide[9].value, xiTide[10].value, xiTide[11].value, xiTide[12].value)
    # println("xiTide: ", xiTide)
    # println("comp: ", comp)
    # println("q: ", q)
    # println("chi1: ", chi1)
    # Compute Kerr BH ISCO radius
    Z1_ISCO = 1.0 + ((1.0 - chi1^2)^(1. /3.))*((1.0+chi1)^(1. /3.) + (1.0-chi1)^(1. /3.))
    Z2_ISCO = sqrt(3.0*chi1^2 + Z1_ISCO^2)
    r_ISCO  = ifelse(chi1>0., 3.0 + Z2_ISCO - sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)), 3.0 + Z2_ISCO + sqrt((3.0 - Z1_ISCO)*(3.0 + Z1_ISCO + 2.0*Z2_ISCO)))
    
    tmpMtorus = alphaTor * xiTide * (1.0-2.0*comp) - betaTor * q*comp * r_ISCO
    
    Mtorus = ifelse(tmpMtorus>0., tmpMtorus, 0.)
    
    
    # Get remnant spin for assumed aligned spin system, from arXiv:1903.11622 Table I and eq. (4), (5) and (6)
    
    p1_remSp = ((-5.44187381e-03*chi1 + 7.91165608e-03) + (2.33362046e-02*chi1 + 2.47764497e-02)*eta)*eta
    p2_remSp = ((-8.56844797e-07*chi1 - 2.81727682e-06) + (6.61290966e-06*chi1 + 4.28979016e-05)*eta)*eta
    p3_remSp = ((-3.04174272e-02*chi1 + 2.54889050e-01) + (1.47549350e-01*chi1 - 4.27905832e-01)*eta)*eta
    
    modelRemSp = (1. + Lambda * p1_remSp + Lambda^2 * p2_remSp) / ((1. + Lambda*p3_remSp^2)*(1. + Lambda*p3_remSp^3))

    modelRemSp = ifelse((chi1 < 0.) & (eta < 0.188), 1., modelRemSp)
    modelRemSp = ifelse(chi1 < -0.5, 1., modelRemSp)
    modelRemSp = ifelse(modelRemSp > 1., 1., modelRemSp)
        
    # Work with spin variables weighted on square of the BH mass over total mass
    S1BH = chi1 * m1ByM^2 
    Shat = S1BH / (m1ByM^2 + m2ByM^2) # this would be = (chi1*m1*m1 + chi2*m2*m2)/(m1*m1 + m2*m2), but chi2=0 by assumption # ANDREA what? chi2 was not put to zero
    
    # Compute fit to L_orb in arXiv:1611.00332 eq. (16)
    #Lorb = (2. *sqrt(3.)*eta + 5.24*3.8326341618708577*eta2 + 1.3*(-9.487364155598392)*eta3)/(1. + 2.88*2.5134875145648374*eta) + ((-0.194)*1.0009563702914628*Shat*(4.409160174224525*eta + 0.5118334706832706*eta2 + (64. - 16. *4.409160174224525 - 4. *0.5118334706832706)*eta3) 
    #+ 0.0851*0.7877509372255369*Shat^2*(8.77367320110712*eta + (-32.060648277652994)*eta2 + (64. - 16. *8.77367320110712 - 4. *(-32.060648277652994))*eta3) + 0.00954*0.6540138407185817*Shat^3*(22.830033250479833*eta + (-153.83722669033995)*eta2 + (64. - 16. *22.830033250479833 - 4. *(-153.83722669033995))*eta3))/(1. + (-0.579)*0.8396665722805308*Shat*(1.8804718791591157 + (-4.770246856212403)*eta + 0. *eta2 + (64. - 64. *1.8804718791591157 - 16. *(-4.770246856212403) - 4. *0.)*eta3)) + 0.3223660562764661*Seta*eta2*(1. + 9.332575956437443*eta)*chi1 + 2.3170397514509933*Shat*Seta*eta3*(1. + (-3.2624649875884852)*eta)*chi1 + (-0.059808322561702126)*eta3*chi12;
    
    Lorb = (2. *sqrt(3.)*eta + 5.24*3.8326341618708577*eta2 + 1.3*(-9.487364155598392)*eta3)/(1. + 2.88*2.5134875145648374*eta) + ((-0.194)*1.0009563702914628*Shat*(4.409160174224525*eta + 0.5118334706832706*eta2 + (64. - 16. *4.409160174224525 - 4. *0.5118334706832706)*eta3) + 0.0851*0.7877509372255369*Shat^2*(8.77367320110712*eta + (-32.060648277652994)*eta2 + (64. - 16. *8.77367320110712 - 4. *(-32.060648277652994))*eta3) + 0.00954*0.6540138407185817*Shat^3*(22.830033250479833*eta + (-153.83722669033995)*eta2 + (64. - 16. *22.830033250479833 - 4. *(-153.83722669033995))*eta3))/(1. + (-0.579)*0.8396665722805308*Shat*(1.8804718791591157 + (-4.770246856212403)*eta + 0. *eta2 + (64. - 64. *1.8804718791591157 - 16. *(-4.770246856212403) - 4. *0.)*eta3)) + 0.3223660562764661*Seta*eta2*(1. + 9.332575956437443*eta)*chi1 + 2.3170397514509933*Shat*Seta*eta3*(1. + (-3.2624649875884852)*eta)*chi1 + (-0.059808322561702126)*eta3*chi12;
    
    #Lorb = (2.*np.sqrt(3.)*eta + 5.24*3.8326341618708577*eta2 + 1.3*(-9.487364155598392)*eta*eta2)/(1. + 2.88*2.5134875145648374*eta) + ((-0.194)*1.0009563702914628*Shat*(4.409160174224525*eta + 0.5118334706832706*eta2 + (64. - 16.*4.409160174224525 - 4.*0.5118334706832706)*eta2*eta) + 0.0851*0.7877509372255369*Shat*Shat*(8.77367320110712*eta + (-32.060648277652994)*eta2 + (64. - 16.*8.77367320110712 - 4.*(-32.060648277652994))*eta2*eta) + 0.00954*0.6540138407185817*Shat*Shat*Shat*(22.830033250479833*eta + (-153.83722669033995)*eta2 + (64. - 16.*22.830033250479833 - 4.*(-153.83722669033995))*eta2*eta))/(1. + (-0.579)*0.8396665722805308*Shat*(1.8804718791591157 + (-4.770246856212403)*eta + 0.*eta2 + (64. - 64.*1.8804718791591157 - 16.*(-4.770246856212403) - 4.*0.)*eta2*eta)) + 0.3223660562764661*Seta*eta2*(1. + 9.332575956437443*eta)*chi1 + 2.3170397514509933*Shat*Seta*eta2*eta*(1. + (-3.2624649875884852)*eta)*chi1 + (-0.059808322561702126)*eta2*eta*chi12;
    
    chif = (Lorb + S1BH)*modelRemSp
    
    # Get remnant mass scaled to a total (initial) mass of 1
    
    p1_remM = ((-1.83417425e-03*chi1 + 2.39226041e-03) + (4.29407902e-03*chi1 + 9.79775571e-03)*eta)*eta
    p2_remM = ((2.33868869e-07*chi1 - 8.28090025e-07) + (-1.64315549e-06*chi1 + 8.08340931e-06)*eta)*eta
    p3_remM = ((-2.00726981e-02*chi1 + 1.31986011e-01) + (6.50754064e-02*chi1 - 1.42749961e-01)*eta)*eta

    modelRemM = (1. + Lambda * p1_remM + Lambda^2 * p2_remM) / ((1. + Lambda*p3_remM^2)*(1. + Lambda*p3_remM^2))
    modelRemM = ifelse((chi1 < 0.) & (eta < 0.188), 1., modelRemM)
    modelRemM = ifelse(chi1 < -0.5, 1., modelRemM)
    modelRemM = ifelse(modelRemM > 1., 1., modelRemM)
    
    # Compute the radiated-energy fit from arXiv:1611.00332 eq. (27)
    #ANDREA why Erad in Phi is different?
    EradNSBH = (((1. + -2.0/3.0*sqrt(2.))*eta + 0.5609904135313374*eta2 + (-0.84667563764404)*eta3 + 3.145145224278187*eta2*eta2)*(1. + 0.346*(-0.2091189048177395)*Shat*(1.8083565298668276 + 15.738082204419655*eta + (16. - 16. *1.8083565298668276 - 4. *15.738082204419655)*eta2) + 0.211*(-0.19709136361080587)*Shat^2*(4.271313308472851 + 0. *eta + (16. - 16. *4.271313308472851 - 4. *0.)*eta2) + 0.128*(-0.1588185739358418)*Shat^3*(31.08987570280556 + (-243.6299258830685)*eta + (16. - 16. *31.08987570280556 - 4. *(-243.6299258830685))*eta2)))/(1. + (-0.212)*2.9852925538232014*Shat*(1.5673498395263061 + (-0.5808669012986468)*eta + (16. - 16. *1.5673498395263061 - 4. *(-0.5808669012986468))*eta2)) + (-0.09803730445895877)*Seta*eta2*(1. + (-3.2283713377939134)*eta)*chi1 + (-0.01978238971523653)*Shat*Seta*eta*(1. + (-4.91667749015812)*eta)*chi1 + 0.01118530335431078*eta3*chi12
    finalMass = (1. -EradNSBH)*modelRemM
    
    # Compute 22 quasi-normal mode dimensionless frequency
    kappaOm = sqrt(log(2. -chif)/log(3.))
    omega_tilde = (1.0 + kappaOm*(1.5578*exp(1im*2.9031) + 1.9510*exp(1im*5.9210)*kappaOm + 2.0997*exp(1im*2.7606)*kappaOm^2 + 1.4109*exp(1im*5.9143)*kappaOm^3 + 0.4106*exp(1im*2.7952)*(kappaOm^4)))
    
    fring = 0.5*real(omega_tilde)/pi/finalMass
    
    rtide = xiTide * (1.0 - 2.0 * comp) / (q*comp)
    
    q_factor = 0.5*real(omega_tilde)/imag(omega_tilde)
    
    ftide = abs(1.0/(pi*(chi1 + sqrt(rtide^3)))*(1.0 + 1.0 / q))
    
    # Now compute last amplitude quantities
    fring_tilde = 0.99 * 0.98 * fring
    
    gamma_correction = ifelse(Lambda > 1.0, 1.25, 1.0 + 0.5*Lambda - 0.25*Lambda^2)
    delta_2_prime = ifelse(Lambda > 1.0, 1.62496*0.25*(1. + tanh(4.0*((ftide/fring_tilde - 1.)-0.0188092)/0.338737)), del2 - 2. *(del2 - 0.81248)*Lambda + (del2 - 0.81248)*Lambda^2)
    
    sigma = delta_2_prime * fring / q_factor
    
    # Determine the type of merger we see and determine coefficients
    epsilon_tide = ifelse(ftide < fring, 0., 2. *0.25*(1 + tanh(4.0*(((ftide/fring_tilde - 1.)*(ftide/fring_tilde - 1.) - 0.571505*comp - 0.00508451*chi1)+0.0796251)/0.0801192)))
    
    epsilon_ins  = ifelse(ftide < fring, ifelse(1.29971 - 1.61724 * (Mtorus + 0.424912*comp + 0.363604*sqrt(eta) - 0.0605591*chi1)>1., 1., 1.29971 - 1.61724 * (Mtorus + 0.424912*comp + 0.363604*sqrt(eta) - 0.0605591*chi1)), ifelse(Mtorus > 0., 1.29971 - 1.61724 * (Mtorus + 0.424912*comp + 0.363604*sqrt(eta) - 0.0605591*chi1), 1.))
    
    sigma_tide   = ifelse(ftide < fring, ifelse(Mtorus>0., 0.137722 - 0.293237*(Mtorus - 0.132754*comp + 0.576669*sqrt(eta) - 0.0603749*chi1 - 0.0601185*chi1^2 - 0.0729134*chi1^3), 0.5*(0.137722 - 0.293237*(Mtorus - 0.132754*comp + 0.576669*sqrt(eta) - 0.0603749*chi1 - 0.0601185*chi1^2 - 0.0729134*chi1^3) + 0.5*(1. - tanh(4.0*(((ftide/fring_tilde - 1.)*(ftide/fring_tilde - 1.) - 0.657424*comp - 0.0259977*chi1)+0.206465)/0.226844)))),0.5*(1. - tanh(4.0*(((ftide/fring_tilde - 1.)*(ftide/fring_tilde - 1.) - 0.657424*comp - 0.0259977*chi1)+0.206465)/0.226844)))
    
    f0_tilde_PN  = ifelse(ftide < fring, ifelse(Mtorus>0., ftide / (M*GMsun_over_c3), ((1.0 - 1.0 / q) * fring_tilde + epsilon_ins * ftide / q)/(M*GMsun_over_c3)), ifelse(Lambda>1., fring_tilde/(M*GMsun_over_c3), ((1.0 - 0.02*Lambda + 0.01*Lambda^2)*0.98*fring)/(M*GMsun_over_c3)))
    
    f0_tilde_PM  = ifelse(ftide < fring, ifelse(Mtorus>0., ftide / (M*GMsun_over_c3), ((1.0 - 1.0 / q) * fring_tilde + ftide/q)/(M*GMsun_over_c3)), ifelse(Lambda>1., fring_tilde/(M*GMsun_over_c3), ((1.0 - 0.02*Lambda + 0.01*Lambda^2)*0.98*fring)/(M*GMsun_over_c3)))
    
    f0_tilde_RD  = ifelse(ftide < fring, 0., ifelse(Lambda>1., fring_tilde/(M*GMsun_over_c3), ((1.0 - 0.02*Lambda + 0.01*Lambda^2)*0.98*fring)/(M*GMsun_over_c3)))

    # This can be used to output the merger type if needed
    
    v = @. (fgrid*pi)^(1. /3.)
    v2 = v .* v
    xdot = @. xdotaN*(v^10)*(1. + xdota2*v2 + xdota3 * fgrid*pi + xdota4 * fgrid*pi*v + xdota5 * v2*fgrid*pi + (xdota6 + xdota6log * 2. *log(v)) * (fgrid*pi)^2 + xdota7 * v*(fgrid*pi)^2)
    ampfacTime = @. sqrt(abs(pi / (1.5 * v * xdot)))
    
    AmpPNre = @. ampfacTime * AN * v2 * (1. + A2*v2 + A3 * fgrid*pi + A4 * v*fgrid*pi + A5 * v2*fgrid*pi + (A6 + A6log * 2. *log(v)) * (fgrid*pi)^2)
    AmpPNim = @. ampfacTime * AN * v2 * (A5imag * v2*fgrid*pi + A6imag * (fgrid*pi)^2)
    
    aPN = @. sqrt(AmpPNre * AmpPNre + AmpPNim * AmpPNim)
    aPM = @. (gamma_correction * g1 * (fgrid^(5. /6.)))

    LRD = @. sigma^2 / ((fgrid - fring) * (fgrid - fring) + sigma^2*0.25)
    aRD = @. epsilon_tide * del1 * LRD * (fgrid^(-7. /6.))
    
    wMinusf0_PN = @. 0.5 * (1. - tanh(4. *(fgrid - (epsilon_ins * f0_tilde_PN)*M*GMsun_over_c3)/(d0 + sigma_tide)))
    wMinusf0_PM = @. 0.5 * (1. - tanh(4. *(fgrid - f0_tilde_PM*M*GMsun_over_c3)/(d0 + sigma_tide)))
    wPlusf0     = @. 0.5 * (1. + tanh(4. *(fgrid - f0_tilde_RD*M*GMsun_over_c3)/(d0 + sigma_tide)))
    
    amplitudeIMR = @. ifelse(fgrid < fcutPar, (aPN * wMinusf0_PN + aPM * wMinusf0_PM + aRD * wPlusf0), 0.)
    
    # Defined as in LALSimulation - LALSimIMRPhenomD.c line 332. Final units are correctly Hz^-1
    Overallamp = 2. * sqrt(5. /(64. *pi)) * M * GMsun_over_c2_Gpc * M * GMsun_over_c3 / dL
    return Overallamp.*amplitudeIMR


end

# deprecated 
# """
# helper function _xiTide_solver(comp, q, chi) for the computation of the parameter :math:`xi_{tide}` in `arXiv:1509.00512 <https://arxiv.org/abs/1509.00512>`_ eq. (8) as a function of the NS compactness, the binary mass ratio and BH spin.
# Used in the PhenomNSBH model.
# """

# function _xiTide_solver(comp, q, chi)

#     # Coefficients of eq. (8) of arXiv:1509.00512, using as variable sqrt(xi) (so order 10 polynomial)
#     mu = q*comp

#     # roots_tmp = roots(Polynomial([-3. *q*mu^2*chi^2, 0., 6. *q*mu, 0., -3. *q, 0., 0., 2. *chi*(mu^(3. /2.)), -3. *mu, 0., 1.]))
#     # roots_tmp_real = @. real(roots_tmp[(abs(imag(roots_tmp))<1e-5) & (real(roots_tmp)>0.)])
#     # res = maximum(roots_tmp_real.^2)
#     res = find_zero( (x,p)-> x^10 - 3*p[2]*x^8 + 2. *p[3]*(p[2]^(3. /2.)) *x^7 - 3. *p[1]*x^4 + 6. *p[1]*p[2]*x^2 - 3. *p[1]*p[2]^2*p[3]^2
#                     ,20., Order1(),[q, mu, chi], maxiters=200)^2
#     return res
# end





end
