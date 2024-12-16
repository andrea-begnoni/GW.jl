#! format: off

#############################################################
# IMRPhenomD WAVEFORM in the TIGER Framework
#############################################################

"""
ToDo: Need documentation
"""
function PolAbs(model::PhenomD_TIGER,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    o1,
    o2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing
)

    #calculate amplitude of waveform
    amp = Ampl(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        o1,
        o2;
        fcutPar = fcutPar,
        fInsJoin_Ampl = fInsJoin_Ampl,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc,
        container = container
    )

    # take into account inclination 
    hp = @. 0.5 * (1.0 + (cos(iota))^2) .* amp
    hc = @. cos(iota) .* amp

    # return plus and cross polarization absolute values
    return [hp, hc]
    
end

"""
ToDo: Need documentation
"""
function Pol(model::PhenomD_TIGER,
    f::AbstractVector,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    iota,
    o1,
    o2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing
)

    hp, hc = PolAbs(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        iota,
        o1,
        o2;
        fcutPar = fcutPar,
        fInsJoin_Ampl = fInsJoin_Ampl,
        GMsun_over_c3 = GMsun_over_c3,
        GMsun_over_c2_Gpc = GMsun_over_c2_Gpc,
        container = container 
    )

    # Return polarization with correct relative phase.
    # Global phase excluded and provided by Phi().
    return [hp, 1im .* hc]

end

"""
For now this is just the PhenomD_NRTidal amplitude
"""
function Ampl(model::PhenomD_TIGER,
    f,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    o1,
    o2;
    fcutPar = 0.2,
    fInsJoin_Ampl = 0.014,
    GMsun_over_c3 = uc.GMsun_over_c3,
    GMsun_over_c2_Gpc = uc.GMsun_over_c2_Gpc,
    container = nothing,
)
    # For now, just use Phenom D
    return Ampl(PhenomD_NRTidal(), f, mc, eta, chi1, chi2, dL, o1, o2, fcutPar = fcutPar,
    fInsJoin_Ampl = fInsJoin_Ampl,
    GMsun_over_c3 = GMsun_over_c3,
    GMsun_over_c2_Gpc = GMsun_over_c2_Gpc)
end

"""
For now this is just the PhenomD_NRtidal case
"""
function Phi(model::PhenomD_TIGER,
    f,
    mc,
    eta,
    chi1,
    chi2,
    o1,
    o2;
    fInsJoin_PHI = 0.018,
    fcutPar = 0.2,
    GMsun_over_c3 = uc.GMsun_over_c3,
    container = nothing,
)
    # For now, just use Phenom D
    return Phi(PhenomD_NRTidal(), f, mc, eta, chi1, chi2, o1, o2, fInsJoin_PHI = fInsJoin_PHI, fcutPar = fcutPar,  GMsun_over_c3 = GMsun_over_c3)
end