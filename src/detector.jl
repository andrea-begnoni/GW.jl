module detector


using ..waveform
import ..UtilsAndConstants as uc     

##### JULIA PACKAGES
using Trapz
using DelimitedFiles
using Interpolations
using ForwardDiff
using LinearAlgebra
using Base.Threads
using BenchmarkTools
using HDF5
using ProgressMeter
using Base.Threads
using Dates


export DetectorStructure, DetectorCoordinates, Detector, _readASD, _readPSD, getCoords, CE1Id_coordinates, CE1Id, CE2NM_coordinates,
         CE2NM, CE2NSW_coordinates, CE2NSW, ETS_coodinates, ETS, ETLS_coodinates, ETLS, ETMR_coordinates, ETMR, ETLMR_coordinates, ETLMR, 
         LIGO_L_coordinates, LIGO_L, LIGO_H_coordinates, LIGO_H, VIRGO_coordinates, VIRGO, KAGRA_coordinates, KAGRA, _available_detectors,
          _define_events, _deltLoc, _patternFunction, PolarizationDet, PhaseDet, Strain, SNR, FisherMatrix, _read_Fishers_SNRs

# Here we define all the structures used inside this module
# define the detector structures
abstract type DetectorStructure end

# this one contains only the coordinates of the detector, the orientation and the arm aperture, all in radians
mutable struct DetectorCoordinates <: DetectorStructure
    latitude_rad::Float64
    longitude_rad::Float64
    orientation_rad::Float64
    arm_aperture_rad::Float64
end

# this one contains the detector structure plus the noise power spectrum and the frequency of the noise, the shape of the detector and a label to identify it when saving the results
mutable struct Detector <: DetectorStructure
    latitude_rad::Float64
    longitude_rad::Float64
    orientation_rad::Float64
    arm_aperture_rad::Float64
    shape::Char
    fNoise::AbstractArray
    psd::AbstractArray
    label::String
end

"""
function to read the amplitude spectral density of the noise
    _readASD("path/to/file.txt")

#### Input arguments:
-  `pathASD` : string, path to the file .txt containing the amplitude spectral density of the noise

#### Optional arguments:
-  `cols` : array, default [1,2], columns of the file containing the frequency and the amplitude spectral density

#### Output:
-  `fNoise` : array, frequency of the noise
-  `psd` : array, amplitude spectral density of the noise

#### Example:
```julia
fNoise, psd = _readASD("path/to/file.txt", [1,2])
```
"""
function _readASD(pathASD::String; cols = [1,2])
    fNoise, psd = open(pathASD) do file

        data = readdlm(file)
        fNoise = data[:, cols[1]]
        psd = data[:, cols[2]] .^ 2

        return fNoise, psd
    end
    return fNoise, psd
end

"""
function to read the power spectral density of the noise, which is the square of the amplitude spectral density
    _readPSD("path/to/file.txt")

#### Input arguments:
-  `pathPSD` : string, path to the file .txt containing the power spectral density of the noise

#### Optional arguments:
-  `cols` : array, default [1,2], columns of the file containing the frequency and the power spectral density

#### Output:
-  `fNoise` : array, frequency of the noise
-  `psd` : array, power spectral density of the noise

#### Example:
```julia
fNoise, psd = _readPSD("path/to/file.txt", [1,2])
```
"""
function _readPSD(pathPSD::String; cols = [1,2])
    fNoise, psd = open(pathPSD) do file

        data = readdlm(file)
        fNoise = data[:, cols[1]]
        psd = data[:, cols[2]]

        return fNoise, psd
    end
    return fNoise, psd
end


# function to extract coordinates of a detector structure
function getCoords(detector)
    return detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad
end

# Here we define the detectors present in the code, to create a new one, just follow the structure of the existing ones, for more information look at the tutorials

# Get the path to the directory of this file
PACKAGE_DIR = @__DIR__

# Go one step back in the path (from ""GW.jl/src" to "GW.jl")
PARENT_DIR = dirname(PACKAGE_DIR)

# Construct the path to the "useful_files" folder from the parent directory
PSDS_DIR = joinpath(PARENT_DIR, "useful_files/")

CE1Id_coordinates=detector.DetectorCoordinates(43.827 * pi/180, -112.825 * pi/180, -45 * pi/180, 90. * pi/180)
CE1Id = detector.Detector(43.827 * pi/180, -112.825 * pi/180, -45 * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"CE_curves/cosmic_explorer.txt")..., "CE1Id")

CE2NM_coordinates=detector.DetectorCoordinates(33.160 * pi/180, -106.480 * pi/180, -105. * pi/180, 90. * pi/180)
CE2NM = detector.Detector(33.160 * pi/180, -106.480 * pi/180, -105. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"CE_curves/cosmic_explorer_20km.txt")..., "CE2NM")

CE2NSW_coordinates=detector.DetectorCoordinates(-34. * pi/180, 145. * pi/180, 0. * pi/180, 90. * pi/180)
CE2NSW=detector.Detector(-34. * pi/180, 145. * pi/180, 0. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"CE_curves/cosmic_explorer_20km.txt")...,  "CE2NSW")

ETS_coodinates = detector.DetectorCoordinates((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 60. * pi/180)
ETS=detector.Detector((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 60. * pi/180, 'T', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETS")

ETLS_coodinates = detector.DetectorCoordinates((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 90. * pi/180)
ETLS = detector.Detector((40. +31. /60.) * pi/180, (9. +25. /60.) * pi/180, 0. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETLS")

ETMR_coordinates = detector.DetectorCoordinates((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 60. * pi/180)
ETMR = detector.Detector((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 60. * pi/180, 'T', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETMR")

ETLMR_coordinates = detector.DetectorCoordinates((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 90. * pi/180)
ETLMR = detector.Detector((50. +43. /60. +23. /3600.) * pi/180, (5. +55. /60. +14. /3600.) * pi/180, 0. * pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"ET_curves/ET-0000A-18_ETDSensitivityCurveTxtFile.txt", cols=[1,4])...,  "ETLMR")

LIGO_L_coordinates = detector.DetectorCoordinates(30.563 * pi/180, -90.774 * pi/180, 242.71636956358617* pi/180, 90. * pi/180)
LIGO_L = detector.Detector(30.563 * pi/180, -90.774 * pi/180, 242.71636956358617* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/aligo_O3actual_L1.txt")...,  "LIGO_L")

LIGO_H_coordinates = detector.DetectorCoordinates(46.455 * pi/180, -119.408 * pi/180, 170.99924234706103* pi/180, 90. * pi/180)
LIGO_H = detector.Detector(46.455 * pi/180, -119.408 * pi/180, 170.99924234706103* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/aligo_O3actual_H1.txt")...,  "LIGO_H")

VIRGO_coordinates = detector.DetectorCoordinates(43.631 * pi/180, 10.504 * pi/180, 115.56756342034298* pi/180, 90. * pi/180)
VIRGO = detector.Detector(43.631 * pi/180, 10.504 * pi/180, 115.56756342034298* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/avirgo_O3actual.txt")..., "VIRGO")

KAGRA_coordinates = detector.DetectorCoordinates(36.412 * pi/180, 137.306 * pi/180, 15.396* pi/180, 90. * pi/180)
KAGRA = detector.Detector(36.412 * pi/180, 137.306 * pi/180, 15.396* pi/180, 90. * pi/180, 'L', _readASD(PSDS_DIR*"observing_scenarios_paper/kagra_128Mpc.txt")...,  "KAGRA")


# with this function you can check all the detectors available, it also prints the length of the arms for future interferomenters
function _available_detectors()
    println("The detectors available are: ")
    return ["CE1Id, 40km", "CE2NM, 20km", "CE2NSW, 20km", "ETS, 10km", "ETLS, 10km", "ETMR, 10km", "ETLMR, 10km", "LIGO_L", "LIGO_H", "VIRGO", "KAGRA"]
end

# this function returns the detector structure given the name of the detector
function _available_detectors(det::String)
    # check if det is among the tabulated detectors
    if det=="CE1Id"
        return CE1Id
    elseif det=="CE2NM"
        return CE2NM
    elseif det=="CE2NSW"
        return CE2NSW
    elseif det=="ETS"
        return ETS
    elseif det=="ETLS"
        return ETLS
    elseif det=="ETMR"
        return ETMR
    elseif det=="ETLMR"
        return ETLMR
    elseif det=="LIGO_L"
        return LIGO_L
    elseif det=="LIGO_H"
        return LIGO_H
    elseif det=="VIRGO"
        return VIRGO
    elseif det=="KAGRA"
        return KAGRA
    else
        println("Detector not available. The ones available are: CE1Id, CE2NM, CE2NSW, ETS (triangular), ETLS, ETMR (triangular), ETLMR, LIGO_L, LIGO_H, VIRGO, KAGRA")
    end
end

# function to check that the event have sense and the variable are ordered correnctly
function _define_events(mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal)

    if mc <= 0 || mc >= 200
        error("Chirp mass out of range, mc = $mc")
    end
    if eta <= 0 || eta > 0.25   
        error("Symmetric mass ratio out of range, eta = $eta")
    end
    if dL <= 0
        error("Luminosity distance out of range, dL = $dL")
    end
    if theta < 0 || theta > pi
        error("Sky position angle out of range, theta = $theta")
    end
    if phi < 0 || phi > 2 * pi
        error("Sky position angle out of range, phi = $phi")
    end
    if iota < 0 || iota > pi
        error("Inclination angle out of range, iota = $iota")
    end
    if psi < 0 || psi > 2 * pi
        error("Polarisation angle out of range, psi = $psi")
    end
    if tcoal < 0
        error("Time of coalescence out of range, tcoal = $tcoal")
    end
    if phiCoal < 0 || phiCoal > 2 * pi
        error("GW phase at coalescence out of range, phiCoal = $phiCoal")
    end
    if chi1 < -1 || chi1 > 1
        error("Spin component out of range, chi1 = $chi1")
    end
    if chi2 < -1 || chi2 > 1
        error("Spin component out of range, chi2 = $chi2")
    end

    println("Everything looks fine")
    return 0
end

"""
The function computes the time (in seconds) it takes for the GW to go from Earth center, which is the point to which all GWs are referred, to the location of the detector.
    _deltLoc(theta, phi, tRef, DetectorCoordinates) 

#### Input arguments:
-  `theta` : float, in rad
-  `phi`   : float, in rad
-  `tRef`  : float, in LMST, it is the time at which the GW arrives at the center of the Earth, it is equal to `tcoal` if we are not considering the Earth motion during the measurement
-  `DetectorCoordinates` : structure, containing the coordinates of the detector

#### Optional arguments:
-  `REarth_km` : float, default 6371.0, Earth radius in km
-  `clight_kms` : float, default 299792.458, speed of light in km/s

#### Output:
- `delt`  : float, in seconds

#### Example:
```julia
delt = _deltLoc(0.1, 0.2, 0.3, CE1Id_coordinates)
```
"""
function _deltLoc(
    theta::Union{Float64,ForwardDiff.Dual},
    phi::Union{Float64,ForwardDiff.Dual},
    tRef,
    DetectorCoordinates::DetectorStructure;
    REarth_km = uc.REarth_km,
    clight_kms = uc.clight_kms,
)

    # Time needed to go from Earth center to detector location

    ra, dec = uc._ra_dec_from_theta_phi_rad(theta, phi)

    comp1 = 
        @. cos(dec) *
        cos(ra) *
        cos(DetectorCoordinates.latitude_rad) *
        cos(DetectorCoordinates.longitude_rad + 2.0 * pi * tRef)
    comp2 =
        @. cos(dec) *
        sin(ra) *
        cos(DetectorCoordinates.latitude_rad) *
        sin(DetectorCoordinates.longitude_rad + 2.0 * pi * tRef)
    comp3 = sin(dec) * sin(DetectorCoordinates.latitude_rad)
    # The minus sign arises from the definition of the unit vector pointing to the source
    delt = @.  -REarth_km * (comp1 + comp2 + comp3) / clight_kms

    return delt # in seconds
end

"""
The function computes the pattern functions of the detector for a GW event, they are function of the detector position at the time of arrival of the GW, the sky position of the source, the polarisation of the GW.

    _patternFunction(theta, phi, psi, tRef, DetectorCoordinates)

#### Input arguments:
-  `theta` : float, in rad
-  `phi`   : float, in rad
-  `psi`   : float, in rad
-  `tRef`  : float, in LMST, it is the time at which the GW arrives at the center of the Earth, it is equal to `tcoal` if we are not considering the Earth motion during the measurement
-  `DetectorCoordinates` : structure, containing the coordinates of the detector

#### Optional arguments:
-  `alpha_grad` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry

#### Output:
- `Fp`  : float, plus pattern function
- `Fc`  : float, cross pattern function

#### Example:
```julia
Fp, Fc = _patternFunction(0.1, 0.2, 0.3, 0.4, CE1Id_coordinates)
```

ToDo: Not up to date
"""
function _patternFunction(
    model::Model,
    DetectorCoordinates::DetectorStructure,
    tRef,
    theta::Union{Float64,ForwardDiff.Dual},
    phi::Union{Float64,ForwardDiff.Dual},
    psi::Union{Float64,ForwardDiff.Dual};
    alpha_grad = 0.0,
)

    alpha_rad = alpha_grad * pi / 180.0
    lambda_rad = DetectorCoordinates.latitude_rad                        # latitude of detector site
    zeta_rad = DetectorCoordinates.arm_aperture_rad                      # angle between intereferrometer arms
    phi_mv_rad = @. DetectorCoordinates.longitude_rad + 2.0 * pi * tRef  # angles between the local meridian and the vernal point
    gamma_rad = DetectorCoordinates.orientation_rad + alpha_rad          # orientation of the detector’s arms with respec to local geographical directions, s measured counterclockwise from East to the bisector of the interferometer arms
    
    # detector tensor components in the celestial sphere frame
    cd = -cos(2*gamma_rad - zeta_rad)/4 + cos(2*gamma_rad + zeta_rad)/4
    sd = -sin(2*gamma_rad - zeta_rad)/4 + sin(2*gamma_rad + zeta_rad)/4   
    dt_cs_1 = @. cd/4*(cos(2*phi_mv_rad)*(3-cos(2*lambda_rad))-cos(2*lambda_rad)-1) - sd*sin(2*phi_mv_rad)*sin(lambda_rad)
    dt_cs_2 = @. cd/4*(3-cos(2*lambda_rad))*sin(2*phi_mv_rad) + sd*sin(lambda_rad)*cos(2*phi_mv_rad)
    dt_cs_3 = @. -cd/2*(sin(2*lambda_rad)*cos(phi_mv_rad)) + sd*sin(phi_mv_rad)*cos(lambda_rad)
    dt_cs_4 = @. cd/4*(cos(2*phi_mv_rad)*(cos(2*lambda_rad)-3)-cos(2*lambda_rad)-1) + sd*sin(2*phi_mv_rad)*sin(lambda_rad)
    dt_cs_5 = @. -cd/2*sin(phi_mv_rad)*sin(2*lambda_rad) - sd*cos(phi_mv_rad)*cos(lambda_rad)

    # calculate transformation matrix from the GW frame to the celestial sphere frame
    ra, dec = uc._ra_dec_from_theta_phi_rad(theta, phi)
    mat_1 = [sin(ra)*cos(psi)-cos(ra)*sin(dec)*sin(psi)   -sin(ra)*sin(psi)-cos(ra)*sin(dec)*cos(psi)   -cos(ra)*cos(dec);
             -sin(ra)*sin(dec)*sin(psi)-cos(ra)*cos(psi)  -sin(ra)*sin(dec)*cos(psi)+sin(psi)*cos(ra)  -sin(ra)*cos(dec);
             sin(psi)*cos(dec)                            cos(dec)*cos(psi)                            -sin(dec)]
    
    # contract polarization tensor with detector tensor in celestial sphere frame
    pol_list = _list_polarizations(model)
    pattern_functions = []
    for pol_name in pol_list
        # poltensor in GW frame
        pol_tensor = uc.polarization_dict[pol_name]
        # transfrom to celestial sphere frame
        pol_tensor = mat_1*pol_tensor*transpose(mat_1)

        # contract detector and polarization tensor

        fp = pol_tensor[1,1] .* dt_cs_1 .+ 2 .*pol_tensor[1,2] .* dt_cs_2 .+ 2 .*pol_tensor[1,3] .* dt_cs_3 .+ pol_tensor[2,2] .* dt_cs_4 .+ 2 .*pol_tensor[2,3] .* dt_cs_5 .- pol_tensor[3,3] .* (dt_cs_1 + dt_cs_4)
             
        push!(pattern_functions, fp)
    end

    return pattern_functions
end

"""
Compute the amplitude of the GW signal projected on the detector tensor for each polarization, given a precomputed list of polarizations, such as outputed by Pol(model, DetectorCoordinates, ... ).

    PolarizationDet(model, DetectorCoordinates, pol, f, mc, eta, dL, theta, phi, iota, psi, tcoal, alpha = 0., useEarthMotion = false)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `pol`: 
    -  `f` : array, frequency of the GW signal, Hz
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `dL` : float, luminosity distance, Gpc
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `iota` : float, inclination angle of the orbital angular momentum to the line of sight toward the detector, radians
    -  `psi` : float, polarisation angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days

    #### Optional arguments:
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular  geometry.
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement

    #### Output:
    - Returns a list of vectors. Each vector corresponds to a contribution of a polarization of the signal, projected onto the detector.

    Example: 
        ToDo: need example
"""
function PolarizationDet(model::Model,
    DetectorCoordinates::DetectorStructure,
    pol::AbstractArray,
    f::AbstractArray,
    mc::Union{Float64,ForwardDiff.Dual},
    eta::Union{Float64,ForwardDiff.Dual},
    dL::Union{Float64,ForwardDiff.Dual},
    theta::Union{Float64,ForwardDiff.Dual},
    phi::Union{Float64,ForwardDiff.Dual},
    iota::Union{Float64,ForwardDiff.Dual},  #Is actually not required in the function. ToDo: Remove
    psi::Union{Float64,ForwardDiff.Dual},
    tcoal::Union{Float64,ForwardDiff.Dual};
    alpha = 0.0,
    useEarthMotion = false
)

    if useEarthMotion   # technique from GWFAST
        tcoalRescaled = tcoal .- waveform._tau_star(model, f, mc, eta) ./ (3600.0 * 24.0)
        tRef =
            tcoalRescaled .+
            _deltLoc(theta, phi, tcoalRescaled, DetectorCoordinates) ./ (3600.0 * 24.0)
    else
        tRef = tcoal + _deltLoc(theta, phi, tcoal, DetectorCoordinates) / (3600.0 * 24.0)
    end

    #Fp, Fc = _patternFunction(theta, phi, psi, tRef, DetectorCoordinates, alpha_grad=alpha)
    #Ap = @. Fp .* pol[1]
    #Ac = @. Fc .* pol[2]

    # project polarizations onto detector
    Fp = _patternFunction(model, DetectorCoordinates, tRef, theta, phi, psi, alpha_grad=alpha)

    n_fp = length(Fp)
    pol_det = Vector{typeof(pol[1])}(undef, n_fp)
    for idx in 1:n_fp
        pol_det[idx] = @. Fp[idx] .* pol[idx] 
    end

    return pol_det

end

"""
Need documentation 
"""
function PolarizationDet(model::Model,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc::Union{Float64,ForwardDiff.Dual},
    eta::Union{Float64,ForwardDiff.Dual},
    chi1::Union{Float64,ForwardDiff.Dual},
    chi2::Union{Float64,ForwardDiff.Dual},
    dL::Union{Float64,ForwardDiff.Dual},
    theta::Union{Float64,ForwardDiff.Dual},
    phi::Union{Float64,ForwardDiff.Dual},
    iota::Union{Float64,ForwardDiff.Dual},
    psi::Union{Float64,ForwardDiff.Dual},
    tcoal::Union{Float64,ForwardDiff.Dual},
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    alpha = 0.0,
    useEarthMotion = false
)

    pol = Pol(
        model,
        f,
        mc,
        eta, 
        chi1,
        chi2,
        dL,
        iota,
        Lambda1,
        Lambda2
    )

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        pol,
        f,
        mc,
        eta,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal;
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    return pol_det

end

"""
ToDo: New documentation for this function
This function computes the phase of the waveform seen by the detector, given a waveform model. It already includes the phase due to the Earth motion.

    PhaseDet(model, DetectorCoordinates, f, mc, eta, chi1, chi2, theta, phi, tcoal, phiCoal, DetectorCoordinates, Lambda1, Lambda2; useEarthMotion = false, phase_precomputation = nothing)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `f` : array, frequency of the GW signal, Hz
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `chi1` : float, dimensionless spin component of the first BH
    -  `chi2` : float, dimensionless spin component of the second BH
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days
    -  `phiCoal` : float, GW phase at coalescence, radians
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `phase_precomputation` : array, default nothing, precomputed phase

    #### Output:
    - `phase`  : array, phase of the GW signal, radians

    #### Example:
    ```julia
    phase = PhaseDet(PhenomD(), CE1Id_coordinates, 1:100, 10.0, 0.25, 0.5, 0.5, 0.1, 0.2, 0.3, 0.4)
    ```

"""
function PhaseDet(
    model::Model,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc::Union{Float64,ForwardDiff.Dual},
    eta::Union{Float64,ForwardDiff.Dual},
    chi1::Union{Float64,ForwardDiff.Dual},
    chi2::Union{Float64,ForwardDiff.Dual},
    theta::Union{Float64,ForwardDiff.Dual},
    phi::Union{Float64,ForwardDiff.Dual},
    tcoal::Union{Float64,ForwardDiff.Dual},
    phiCoal::Union{Float64,ForwardDiff.Dual},
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    useEarthMotion = false,
    phase_precomputation = nothing,
)

    # Phase of the GW signal

    # Full GW strain expression (complex)
    if useEarthMotion # technique from GWFAST
        tcoalRescaled = tcoal .- waveform._tau_star(model, f, mc, eta) ./ (3600.0 * 24.0)  #tcoal is in fraction of days, tau_star in seconds

    else
        tcoalRescaled = tcoal
    end

    phiD = (2.0 * pi .* f) .* _deltLoc(theta, phi, tcoalRescaled, DetectorCoordinates) # phase due to Earth motion

    return @. 2.0 * pi * (tcoal * 3600.0 * 24.0) .* f .- phiCoal  .+ phiD #.- Phi


end


"""
ToDo: Documentation
"""
function Strain(model::Model,
    DetectorCoordinates::DetectorStructure,
    pol::AbstractArray,
    phase::AbstractArray,
    f::AbstractArray,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    phiCoal,
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    useEarthMotion = false
)

    phase_det = PhaseDet(
        model,
        DetectorCoordinates,
        f,
        mc,
        eta,
        chi1,
        chi2,
        theta,
        phi,
        tcoal,
        phiCoal,
        Lambda1,
        Lambda2,
        useEarthMotion = useEarthMotion
    )
        
    return sum(pol) .* exp.(1im .* (phase_det .- phase))

end

"""
ToDo: Documentation
This function computes the full strain (complex) as a function of the parameters, at given frequencies, at detector location, as measured by the detector, since it include the pattern functions, via PolarizationDet and PhaseDet.

    Strain(model, DetectorCoordinates,  f, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1, Lambda2, useEarthMotion=false, alpha=0.)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `f` : array, frequency of the GW signal, Hz
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `chi1` : float, dimensionless spin component of the first BH
    -  `chi2` : float, dimensionless spin component of the second BH
    -  `dL` : float, luminosity distance, Gpc
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `iota` : float, inclination angle of the orbital angular momentum to the line of sight toward the detector, radians
    -  `psi` : float, polarisation angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days
    -  `phiCoal` : float, GW phase at coalescence, radians
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry
    -  `ampl_precomputation` : array, default nothing, precomputed amplitude
    -  `phase_precomputation` : array, default nothing, precomputed phase

    #### Output:
    - `strain`  : array, complex strain as seen by the detector

    #### Example:
    ```julia
    strain = Strain(PhenomD(), CE1Id_coordinates, 1:100, 10.0, 0.25, 0.5, 0.5, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5)
    ```
"""
function Strain(model::Model,
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    phiCoal,
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    useEarthMotion = false,
    alpha = 0.0
)

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        Lambda1,
        Lambda2,
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    phase_wave = Phi(
        model,
        f,
        mc,
        eta,
        chi1,
        chi2,
        Lambda1,
        Lambda2,
    )

    strain_det = Strain(
        model,
        DetectorCoordinates,
        pol_det,
        phase_wave,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        phiCoal,
        Lambda1,
        Lambda2,
        useEarthMotion = useEarthMotion
    )
        
    return strain_det

end

"""
ToDo: Need documentatiation
"""
function Strain(model::Union{PhenomHM,PhenomXHM},
    DetectorCoordinates::DetectorStructure,
    f::AbstractArray,
    mc,
    eta,
    chi1,
    chi2,
    dL,
    theta,
    phi,
    iota,
    psi,
    tcoal,
    phiCoal,
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    useEarthMotion = false,
    alpha = 0.0
)

    pol_det = PolarizationDet(
        model,
        DetectorCoordinates,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        Lambda1,
        Lambda2,
        alpha = alpha,
        useEarthMotion = useEarthMotion
    )

    # ToDo: Not happy with this at all
    phase_wave = f * 0.0

    strain_det = Strain(
        model,
        DetectorCoordinates,
        pol_det,
        phase_wave,
        f,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        phiCoal,
        Lambda1,
        Lambda2,
        useEarthMotion = useEarthMotion
    )
        
    return strain_det

end


"""
Compute the *signal-to-noise-ratio*, SNR, as a function of the parameters of the event, as measured a single detector. 
The SNR is computed as the square root of the integral of the signal-to-noise ratio squared.
The integral is computed using the trapezoidal rule.
This function is the main function to compute the SNR, then using multiple-dispatch this function is called when considering a network of detectors and more than a single event.

SNR(model, detector, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, Lambda1=0.0, Lambda2=0.0; fmin=2.0, fmax=nothing, res=1000, useEarthMotion=false, ampl_precomputation=nothing)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `DetectorCoordinates` : structure, containing the coordinates of the detector
    -  `f` : array, frequency of the GW signal, Hz
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `chi1` : float, dimensionless spin component of the first BH
    -  `chi2` : float, dimensionless spin component of the second BH
    -  `dL` : float, luminosity distance, Gpc
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `iota` : float, inclination angle of the orbital angular momentum to the line of sight toward the detector, radians
    -  `psi` : float, polarisation angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `fmin` : float, default 2.0, minimum frequency
    -  `fmax` : float, default nothing, maximum frequency, otherwise the code takes fcut (from _fcut) as fmax
    -  `res` : int, default 1000, resolution of the frequency grid
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `ampl_precomputation` : array, default nothing, precomputed amplitude

    #### Output:
    - `SNR`  : float, signal-to-noise ratio

    #### Example:
    ```julia
    SNR = SNR(PhenomD(), CE1Id , 10.0, 0.25, 0.5, 0.5, 1.0, 0.1, 0.2, 0.3, 0.4)
    ```
"""
function SNR(model::Model,
    detector::Detector,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    Lambda1=0.0,
    Lambda2=0.0;
    fmin=2.0,
    fmax = nothing,
    res = 1000,
    useEarthMotion = false,
    ampl_precomputation = nothing,
)

    
    if isnothing(fmax)
        fcut = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
    else
        fcut_tmp = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
        fcut = ifelse(fcut_tmp > fmax, fmax, fcut_tmp)
    end


    fgrid = 10 .^ (range(log10(fmin), log10(fcut), length = res))
    # Out of the PSD range, we use a constant value of 1, which results in completely negligible contributions
    psdGrid = linear_interpolation(detector.fNoise, detector.psd, extrapolation_bc = 1.0)(fgrid)
    detectorCoordinates = DetectorCoordinates(
        detector.latitude_rad,
        detector.longitude_rad,
        detector.orientation_rad,
        detector.arm_aperture_rad
    )

        if ampl_precomputation === nothing
           
            ampl_precomputation = PolAbs(
                model,
                f,
                mc,
                eta,
                chi1,
                chi2,
                dL,
                iota,
                Lambda1,
                Lambda2 
            )
        else
            ampl_precomputation = ampl_precomputation
        end

    if detector.shape == 'L' # in Julia '' indicates a char, while "" a string 
        Aps, Acs = PolarizationDet(
            model,
            detectorCoordinates,
            ampl_precomputation,
            fgrid,
            mc,
            eta,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            useEarthMotion = useEarthMotion,
        )
        Atot = Aps .* Aps .+ Acs .* Acs
        SNR = 2.0 * sqrt(trapz(fgrid, Atot ./ psdGrid))

    elseif detector.shape == 'T'
        SNRsq_Tshape = zeros(3)

        
        # The signal in 3 arms sums to zero for geometrical reasons, so we can use this to skip some calculations

        Aps1, Acs1 = PolarizationDet(
            model,
            detectorCoordinates,
            ampl_precomputation,
            fgrid,
            mc,
            eta,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            alpha = 0.0,
            useEarthMotion = useEarthMotion,
        )
        Atot1 = Aps1 .* Aps1 .+ Acs1 .* Acs1
        Aps2, Acs2 = PolarizationDet(
            model,
            detectorCoordinates,
            ampl_precomputation,
            fgrid,
            mc,
            eta,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            alpha = 60.0,
            useEarthMotion = useEarthMotion,
        )
        Atot2 = Aps2 .* Aps2 .+ Acs2 .* Acs2
        Aps3, Acs3 = -(Aps1 .+ Aps2), -(Acs1 .+ Acs2)
        Atot3 = Aps3 .* Aps3 .+ Acs3 .* Acs3
        SNRsq_Tshape[1] = trapz(fgrid, Atot1 ./ psdGrid)
        SNRsq_Tshape[2] = trapz(fgrid, Atot2 ./ psdGrid)
        SNRsq_Tshape[3] = trapz(fgrid, Atot3 ./ psdGrid)


        SNR = 2.0 * sqrt(sum(SNRsq_Tshape)) # The factor of two arises by cutting the integral from 0 to infinity
    
    end

    return SNR

end

"""
This function computes the *signal-to-noise-ratio*, SNR, as a function of the parameters of the event, as measured by a NETWORK of detectors.
It relies on the function SNR(..., detector::Detector, ...), which computes the SNR for a single detector.
where the dots indicate the parameters equal to the previous function call.

#### Optional arguments:
-  `precomputation` : bool, default true, if true the waveform is called just once and then reused for each detector.

#### Example:
```julia
    SNR(... [CE1Id, CE2NM], ...)
```

"""
function SNR(model::Model,
    detector::Vector{Detector},
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    Lambda1 = 0.0,
    Lambda2 = 0.0;
    fmin=2.0,
    fmax = nothing,
    res = 1000,
    useEarthMotion = false,
    precomputation = true,
)
##########################
## This part is to precompute the amplitude of the waveform which is the longest part of the computation
    if isnothing(fmax)
        fcut = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
    else
        fcut_tmp = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
        fcut = ifelse(fcut_tmp > fmax, fmax, fcut_tmp)
    end


    fgrid = 10 .^ (range(log10(fmin), log10(fcut), length = res))

    if precomputation == true 
        ampl_precomputation = waveform.PolAbs(
            model,
            fgrid,
            mc,
            eta,
            chi1,
            chi2,
            dL,
            iota,
            Lambda1,
            Lambda2;
        )
    else 
        ampl_precomputation = nothing
    end

########################
    SNRList = Vector{Float64}(undef,length(detector))
    for i in 1:length(detector)
        SNRList[i] = SNR(
            model,
            detector[i],
            mc,
            eta,
            chi1,
            chi2,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            Lambda1,
            Lambda2,
            fmin=fmin,
            fmax=fmax,
            res=res,
            useEarthMotion = useEarthMotion,
            ampl_precomputation = ampl_precomputation
        )
    end

    return sqrt(sum(SNRList.^2))
end

"""
This function computes the *signal-to-noise-ratio*, SNR, as a function of the parameters of an ARRAY of event, as measured by a SINGLE detector or a NETWORK of detectors.
It relies on the function SNR(..., detector::Detector, ...), which computes the SNR for a single detector.
where the dots indicate the parameters equal to the previous function call.

It is possible to save the SNRs in a file, if the optional argument `auto_save` is set to true. Use `name_folder` : string, to decide the name of the folder where the results are saved,
if the default is left it saves BBH in the folder "output/BBH" and so on for each source type. The file is saved in the folder `output/name_folder/SNRs.h5` 
and contains the SNRs for the events in the catalog. It contains also the parameters of the events if the optional argument `save_catalog` is set to true.
"""
function SNR(model::Model,
    detector::Union{Detector, Vector{Detector}},
    mc::AbstractArray,
    eta::AbstractArray,
    chi1::AbstractArray,
    chi2::AbstractArray,
    dL::AbstractArray,
    theta::AbstractArray,
    phi::AbstractArray,
    iota::AbstractArray,
    psi::AbstractArray,
    tcoal::AbstractArray,
    Lambda1=0.0,
    Lambda2=0.0;
    fmin=2.0,
    fmax = nothing,
    res = 1000,
    auto_save = false,
    name_folder = "BBH",
    save_catalog = false,
    useEarthMotion = false,
    precomputation = true,
)
    nEvents = length(mc)
    SNRs = Vector{Float64}(undef, nEvents)

    if typeof(model) == PhenomD || typeof(model) == PhenomHM || typeof(model) == PhenomXAS || typeof(model) == PhenomXHM
        Lambda1 = zeros(nEvents)
        Lambda2 = zeros(nEvents)

    elseif typeof(model) == PhenomNSBH
        Lambda2 = zeros(nEvents)
        if name_folder == "BBH"
            name_folder = "NSBH"
        end

    elseif typeof(model) == TaylorF2
        if Lambda1 == 0.
            Lambda1 = zeros(nEvents)
        end
        if Lambda2 == 0.
            Lambda2 = zeros(nEvents)
        end
        
    else
        if name_folder == "BBH"
            name_folder = "BNS"
        end

    end

    elapsed_time = @elapsed @showprogress desc="Computing SNRs..."  @threads for ii in 1:nEvents  
                    SNRs[ii]=SNR(
                        model,
                        detector,
                        mc[ii], 
                        eta[ii], 
                        chi1[ii], 
                        chi2[ii], 
                        dL[ii], 
                        theta[ii], 
                        phi[ii], 
                        iota[ii], 
                        psi[ii], 
                        tcoal[ii],
                        Lambda1[ii],
                        Lambda2[ii], 
                        fmin=fmin, 
                        fmax=fmax, 
                        res = res, 
                        useEarthMotion = useEarthMotion,
                        precomputation = precomputation)
    end 

    println("SNRs computed!")
    if elapsed_time > 60.0
        elapsed_time = elapsed_time/60
        println("The evaluation took: ", elapsed_time, " minutes.")
    else
        println("The evaluation took: ", elapsed_time, " seconds.")
    end

    if auto_save == true
        path = pwd()
        mkpath("output/"*name_folder)
        path = pwd()*"/output/"*name_folder*"/"
        date = Dates.now()
        date_format = string(Dates.format(date, "e dd u yyyy HH:MM:SS"))
        h5open(path*"SNRs.h5", "w") do file
            attributes(file)["number_events"] = nEvents
            if typeof(detector) == Vector{Detector}
                label = [detector[i].label for i in eachindex(detector)]
                attributes(file)["Detectors"] = label
                attributes(file)["What_this_file_contains"] = "This file contains the SNR for the "*string(nEvents)*" events of the catalog 
                                generated with the "*string(typeof(model))*" waveform model. The SNR is calculated for the "*join(label,",")*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                attributes(file)["date"] = date_format
            else
                attributes(file)["Detectors"] = detector.label
                attributes(file)["What_this_file_contains"] = "This file contains the SNR for the "*string(nEvents)*" events of the catalog 
                                generated with the "*string(typeof(model))*" waveform model. The SNR is calculated for the "*detector.label*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                attributes(file)["date"] = date_format
            end   
            
            write(file, "SNRs", SNRs) 
            if save_catalog 
                write(file, "mc", mc)
                write(file, "eta", eta)
                write(file, "chi1", chi1)
                write(file, "chi2", chi2)
                write(file, "dL", dL)
                write(file, "theta", theta)
                write(file, "phi", phi)
                write(file, "iota", iota)
                write(file, "psi", psi)
                write(file, "tcoal", tcoal)
                write(file, "Lambda1", Lambda1)
                write(file, "Lambda2", Lambda2)
            end
        end
    end

    return SNRs
end

"""
This function computes the Fisher Matrix for a single detector, given the parameters of the event and the detector.
To do the computation it uses the function FisherMatrix_internal(...) for L-shaped detectors and FisherMatrix_Tdetector(...) for T-shaped detectors.
Thus it is a wrapper function that calls the correct function depending on the shape of the detector. More information on the Fisher Matrix computation can be found in the documentation of FisherMatrix_internal(...) and FisherMatrix_Tdetector(...).

    FisherMatrix(model, detector , mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1=0.0, Lambda2=0.0, res=1000, useEarthMotion=false, alpha=0.0, rho_thres=12., fmin=2., fmax=nothing, coordinate_shift=true, return_SNR=false)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `detector` : structure, containing the detector information
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `chi1` : float, dimensionless spin component of the first BH
    -  `chi2` : float, dimensionless spin component of the second BH
    -  `dL` : float, luminosity distance, Gpc
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `iota` : float, inclination angle of the orbital angular momentum to the line of sight toward the detector, radians
    -  `psi` : float, polarisation angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days
    -  `phiCoal` : float, GW phase at coalescence, radians
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `res` : int, default 1000, resolution of the frequency grid
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry
    -  `rho_thres` : float, default 12., SNR threshold for the computation of the Fisher Matrix
    -  `fmin` : float, default 2.0, minimum frequency
    -  `fmax` : float, default nothing, maximum frequency, otherwise the code takes fcut (from _fcut) as fmax
    -  `coordinate_shift` : bool, default true, valid for T detectors, if true the codes shifts the coordinates of the detector from the center of the triangle to the center of the arms (more realistic scenario, recommended)
    -  `return_SNR` : bool, default false, if true the function returns the SNR of the event (skipping the need to call the SNR function)

    #### Output:
    - `FisherMatrix`  : matrix, Fisher Matrix

    #### Example:
    ```julia
    FisherMatrix = FisherMatrix(PhenomD(), CE1Id , 10.0, 0.25, 0.5, 0.5, 1.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6)
    ```


"""
function FisherMatrix(model::Model,
    detector::Detector,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    Lambda1=0.,
    Lambda2=0.;
    res = 1000,
    useEarthMotion = false,
    alpha = 0.0,
    rho_thres = 12.,
    fmin=2.,
    fmax=nothing,
    coordinate_shift = true,
    return_SNR = false,
)
    #function that is used only to divide between L and T detectors

    if detector.shape =='L'

        return FisherMatrix_internal(
            model,
            detector,
            mc,
            eta,
            chi1,
            chi2,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            phiCoal,
            Lambda1,
            Lambda2,
            rho_thres=rho_thres,
            res = res,
            useEarthMotion = useEarthMotion,
            alpha = alpha,
            fmin=fmin,
            fmax=fmax,
            return_SNR = return_SNR,
        )
    elseif detector.shape =='T'

        return FisherMatrix_Tdetector(
            model,
            detector,
            mc,
            eta,
            chi1,
            chi2,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            phiCoal,
            Lambda1,
            Lambda2,
            rho_thres=rho_thres,
            res = res,
            useEarthMotion = useEarthMotion,
            alpha = alpha,
            fmin=fmin,
            fmax=fmax,
            coordinate_shift=coordinate_shift,
            return_SNR = return_SNR,
        )
    else
        error("detector shape not recognized")
    end
end

"""
This is the function that computes the Fisher Matrix for a single detector with L shape.
It computes the derivatives of the strain w.r.t. each parameter and then computes the Fisher Matrix. It is called by the function FisherMatrix and FisherMatrix_Tdetector.
There is no need to call this function in your computations since all the logic is implemented in the FisherMatrix function. 

"""

function FisherMatrix_internal(model::Model,
    detector::Detector,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    Lambda1=0.,
    Lambda2=0.;
    res = 1000,
    useEarthMotion = false,
    alpha = 0.0,
    rho_thres = 12.,
    fmin=2.,
    fmax=nothing,
    return_SNR = false,
)

    if model == TaylorF2()
        nPar = _npar(model, Lambda1, Lambda2)
    else
        nPar = _npar(model)
    end

    if isnothing(fmax)
        fcut = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
    else
        fcut_tmp = waveform._fcut(model, mc, eta, Lambda1, Lambda2)
        fcut = ifelse(fcut_tmp > fmax, fmax, fcut_tmp)
    end


    fgrid = 10 .^ (range(log10(fmin), log10(fcut), length = res))
    psdGrid = linear_interpolation(detector.fNoise, detector.psd, extrapolation_bc = 1.0)(fgrid)  
    
    # compute SNR and procede only if it is above the threshold
    if rho_thres !==nothing
        SNRval = SNR(
            model,
            detector,
            mc,
            eta,
            chi1,
            chi2,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            Lambda1,
            Lambda2,
            fmin = fmin,
            fmax = fmax,
            res = res,
            #ampl_precomputation = ampl_precomputation,
        )
        if SNRval < rho_thres
            if return_SNR == true
                return zeros(nPar, nPar), SNRval
            else
                return zeros(nPar, nPar)
            end
        end
    end

    detectorCoordinates = DetectorCoordinates(
        detector.latitude_rad,
        detector.longitude_rad,
        detector.orientation_rad,
        detector.arm_aperture_rad
    )
    
    ###########  Derivatives of the strain w.r.t. each parameter
    strainAutoDiff_real = Matrix{Float64}(undef, res, nPar)
    strainAutoDiff_imag = Matrix{Float64}(undef, res, nPar)
    
    event_parameter = [mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal, Lambda1, Lambda2]
    event_parameter = event_parameter[1:nPar] # cut off non-required parameter,  
    strainAutoDiff_real = ForwardDiff.jacobian(
        x -> real(
            Strain(
                model,
                detectorCoordinates,
                fgrid,
                x... ,
                alpha = alpha,
                useEarthMotion = useEarthMotion
            ),
        ),
        event_parameter,
    )
    strainAutoDiff_imag = ForwardDiff.jacobian(
        x -> imag(
            Strain(
                model,
                detectorCoordinates,
                fgrid,
                x... ,
                alpha = alpha,
                useEarthMotion = useEarthMotion,
            ),
        ),
        event_parameter,
    )

    # It can happen that a certain frequency gives a Nan value, in this case we set the derivative to zero,
    # this happens less than one time per event and usually at the end of the frequency grid.
    strainAutoDiff_real[isnan.(strainAutoDiff_real)] .= 0.0
    strainAutoDiff_imag[isnan.(strainAutoDiff_imag)] .= 0.0
    ######### end of derivatives
    jacobian = Matrix{ComplexF64}(undef, nPar, res)
    for ii in 1:nPar
        jacobian[ii, :] = strainAutoDiff_real[:, ii] + 1im * strainAutoDiff_imag[:, ii]
    end
    jacobian[10,:] /= (3600.0 * 24.0)   # Change the units of the tcoal derivative from days to seconds (this improves conditioning)


    # compute the Fisher matrix
    Fisher = Matrix{Float64}(undef, nPar, nPar)

    for alpha = 1:nPar
        for beta = alpha:nPar
            Fisher[alpha, beta] =
                4.0 *
                trapz(fgrid, real(jacobian[alpha, :] .* conj(jacobian[beta, :])) ./ psdGrid)
            Fisher[beta, alpha] = Fisher[alpha, beta]
        end
    end

    if return_SNR == true
        return Fisher, SNRval
    else
        return Fisher
    end

end

"""
This function computes the *Fisher Matrix*, as a function of the parameters of the event, as measured by a NETWORK of detectors.
It relies on the function FisherMatrix(..., detector::Detector, ...), which computes the Fisher for a single detector.
where the dots indicate the parameters equal to the previous function call.

#### Example:
```julia
    FisherMatrix(... [CE1Id, CE2NM], ...)
```

"""
function FisherMatrix(model::Model,
    detector::Vector{Detector},
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    Lambda1=0.,
    Lambda2=0.;
    res = 1000,
    useEarthMotion = false,
    rho_thres=12.,
    alpha = 0.0,
    fmin=2.0,
    fmax = nothing,
    coordinate_shift = true,
    return_SNR = false,
)

    # compute SNR and procede only if it is above the threshold

    if model == TaylorF2()
        nPar = _npar(model, Lambda1, Lambda2)
    else
        nPar = _npar(model)
    end

    if rho_thres !==nothing
        SNRval = SNR(
            model,
            detector,
            mc,
            eta,
            chi1,
            chi2,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            Lambda1,
            Lambda2,
            fmin = fmin,
            fmax = fmax,
            res = res,
            useEarthMotion = useEarthMotion,
        )
        if SNRval < rho_thres
            if return_SNR == true
                return zeros(nPar, nPar), SNRval
            else
                return zeros(nPar, nPar)
            end
        
        end
    end

    fisherList = Vector{Matrix{Float64}}(undef, length(detector))
    for i in eachindex(detector)
        

        if detector[i].shape == 'L'
            F = FisherMatrix_internal(
                model,
                detector[i],
                mc,
                eta,
                chi1,
                chi2,
                dL,
                theta,
                phi,
                iota,
                psi,
                tcoal,
                phiCoal,
                Lambda1,
                Lambda2,
                rho_thres=nothing,
                res = res,
                useEarthMotion = useEarthMotion,
                alpha = alpha,
                fmin=fmin,
                fmax=fmax,
                return_SNR=false,
            )
        elseif detector[i].shape == 'T'
            F = FisherMatrix_Tdetector(
                model,
                detector[i],
                mc,
                eta,
                chi1,
                chi2,
                dL,
                theta,
                phi,
                iota,
                psi,
                tcoal,
                phiCoal,
                Lambda1,
                Lambda2,
                rho_thres=nothing,
                res = res,
                useEarthMotion = useEarthMotion,
                alpha = alpha,
                fmin=fmin,
                fmax=fmax,
                coordinate_shift = coordinate_shift,
                return_SNR=false,
            )
        end
        fisherList[i] = F

    end
    if return_SNR == true
        return sum(fisherList, dims = 1)[1], SNRval
    else
        return sum(fisherList, dims = 1)[1]
    end

end

"""
This is a helper function that computes the Fisher Matrix for a single detector with T shape. It is called by the function FisherMatrix and calls FisherMatrix_internal.
"""
function FisherMatrix_Tdetector(model::Model,
    detector::Detector,
    mc::Float64,
    eta::Float64,
    chi1::Float64,
    chi2::Float64,
    dL::Float64,
    theta::Float64,
    phi::Float64,
    iota::Float64,
    psi::Float64,
    tcoal::Float64,
    phiCoal::Float64,
    Lambda1=0.,
    Lambda2=0.;
    res = 1000,
    useEarthMotion = false,
    rho_thres=12.,
    alpha = 0.0,
    fmin=2.0,
    fmax = nothing,
    REarth_km = uc.REarth_km,
    coordinate_shift = true,
    return_SNR = false,    
)

    if rho_thres !==nothing

        if model == TaylorF2()
            nPar = _npar(model, Lambda1, Lambda2)
        else
            nPar = _npar(model)
        end

        SNRval = SNR(
            model,
            detector,
            mc,
            eta,
            chi1,
            chi2,
            dL,
            theta,
            phi,
            iota,
            psi,
            tcoal,
            Lambda1,
            Lambda2,
            fmin = fmin,
            fmax = fmax,
            res = res,
            #ampl_precomputation = ampl_precomputation
        )
        if SNRval < rho_thres
            if return_SNR == true
                return zeros(nPar, nPar), SNRval
            else
                return zeros(nPar, nPar)
            end
        end
    end



        # We write the T-detector as three detectors in slightly different positions
        # Detector has to be in the plane orthogonal to the radius
        # write coordinates on a tangent plane (https://en.wikipedia.org/wiki/Local_tangent_plane_coordinates)
        # Note the minus sign because theta goes from north to south

    if coordinate_shift == true
        lat = detector.latitude_rad
        long = detector.longitude_rad
        first = [- cos(lat)*cos(long), 
                - cos(lat)*sin(long), 
                sin(lat)]
        #No sign needed here because phi is counterclosk wise (so goes towards east)
        second = [-sin(long), cos(long), 0.0]
        ETarm           = 10e3                      ## ET arms in m
        ET_cartesian   =      REarth_km * 1e3 .* [sin(lat)*cos(long), sin(lat)*sin(long), cos(lat)]  ## ET center in m
        #These are the positions of the 3 detectors
        ET1_cartesian = @. ET_cartesian + ETarm * (-.5 * first - 0.28867513 * second)
        ET2_cartesian = @. ET_cartesian + ETarm * (+.5 * first - 0.28867513 * second)
        ET3_cartesian = @. ET_cartesian + ETarm * (+0.57735027 * second)

        # go back to spherical coordinates and discard radius (it is REarth at 1e-7)

        ET1_coo = [acos(ET1_cartesian[3]/sqrt(sum(ET1_cartesian.^2))), atan(ET1_cartesian[2], ET1_cartesian[1])]
        ET2_coo = [acos(ET2_cartesian[3]/sqrt(sum(ET2_cartesian.^2))), atan(ET2_cartesian[2], ET2_cartesian[1])]
        ET3_coo = [acos(ET3_cartesian[3]/sqrt(sum(ET3_cartesian.^2))), atan(ET3_cartesian[2], ET3_cartesian[1])]

        ET1 = Detector(ET1_coo[1], ET1_coo[2], detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET2 = Detector(ET2_coo[1], ET2_coo[2], detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET3 = Detector(ET3_coo[1], ET3_coo[2], detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
    else
        ET1 = Detector(detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET2 = Detector(detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
        ET3 = Detector(detector.latitude_rad, detector.longitude_rad, detector.orientation_rad, detector.arm_aperture_rad, 'L', detector.fNoise, detector.psd, detector.label)
    end

    F1 = FisherMatrix_internal(
        model,
        ET1,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        phiCoal,
        Lambda1,
        Lambda2,
        res = res,
        useEarthMotion = useEarthMotion,
        rho_thres = nothing,
        alpha = 0.0,
        fmin=fmin,
        fmax=fmax,
        return_SNR=false,
    )
    F2 = FisherMatrix_internal(
        model,
        ET2,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        phiCoal,
        Lambda1,
        Lambda2,
        res = res,
        useEarthMotion = useEarthMotion,
        rho_thres = nothing,
        alpha = 60.0,
        fmin=fmin,
        fmax=fmax,
        return_SNR=false,
    )
    F3 = FisherMatrix_internal(
        model,
        ET3,
        mc,
        eta,
        chi1,
        chi2,
        dL,
        theta,
        phi,
        iota,
        psi,
        tcoal,
        phiCoal,
        Lambda1,
        Lambda2,
        res = res,
        useEarthMotion = useEarthMotion,
        rho_thres = nothing,
        alpha = 120.0,
        fmin=fmin,
        fmax=fmax,
        return_SNR=false,

    )
    if return_SNR == true
        return F1 + F2 + F3, SNRval
    else
        return F1 + F2 + F3
    end
end

"""
The main function of the code, it computes the Fisher Matrix for an array of events, as measured by a single detector or a network of detectors. It also calculates the 
SNRs if requested and saves the results (Fisher matrices and SNRs) in a file if the optional argument `auto_save` is set to true. The file is saved in the folder `output/name_folder/Fishers_SNRs.h5`


    FisherMatrix(model, detector, mc, eta, chi1, chi2, dL, theta, phi, iota, psi, tcoal, phiCoal,  Lambda1=0.0, Lambda2=0.0, res=1000, useEarthMotion=false, alpha=0.0, rho_thres=12., fmin=2., fmax=nothing, coordinate_shift=true, return_SNR=false)

    #### Input arguments:
    -  `model` : structure, containing the waveform model
    -  `detector` : structure, containing the detector information
    -  `mc` : float, chirp mass, solar masses
    -  `eta` : float, symmetric mass ratio
    -  `chi1` : float, dimensionless spin component of the first BH
    -  `chi2` : float, dimensionless spin component of the second BH
    -  `dL` : float, luminosity distance, Gpc
    -  `theta` : float, sky position angle, radians
    -  `phi` : float, sky position angle, radians
    -  `iota` : float, inclination angle of the orbital angular momentum to the line of sight toward the detector, radians
    -  `psi` : float, polarisation angle, radians
    -  `tcoal` : float, time of coalescence, GMST, fraction of days
    -  `phiCoal` : float, GW phase at coalescence, radians
    -  `Lambda1` : float, tidal parameter of the first object, default 0.0
    -  `Lambda2` : float, tidal parameter of the second object, default 0.0

    #### Optional arguments:
    -  `res` : int, default 1000, resolution of the frequency grid
    -  `useEarthMotion` : bool, default false, if true the Earth motion is considered during the measurement
    -  `alpha` : float, default 0.0, further rotation of the interferometer with respect to the east-west direction, needed for the triangular geometry
    -  `rho_thres` : float, default 12., SNR threshold for the computation of the Fisher Matrix
    -  `fmin` : float, default 2.0, minimum frequency
    -  `fmax` : float, default nothing, maximum frequency, otherwise the code takes fcut (from _fcut) as fmax
    -  `coordinate_shift` : bool, default true, valid for T detectors, if true the codes shifts the coordinates of the detector from the center of the triangle to the center of the arms (more realistic scenario, recommended)
    -  `return_SNR` : bool, default false, if true the function returns the SNR of the event (skipping the need to call the SNR function)
    -  `auto_save` : bool, default false, if true the function saves the results in a file
    -  `name_folder` : string, name of the folder where the results are saved, if the default is left, it saves BBH in the folder "output/BBH" and so on for each source type

    #### Output:
    - `FisherMatrix`  : matrix, Fisher Matrix

    #### Example:
    ```julia
    FisherMatrix = FisherMatrix(PhenomD(), [10.0, 20.], [0.25, 0.25], [0.5, 1.], [0.5, -1.], [1.0, 3.], [0.1, 0.2], [0.2, 0.3], [0.3, 0.4], [0.4, 0.5], [0.5, 0.6], [0.6, 0.7], CE1Id)
    ```


"""
function FisherMatrix(model::Model,
    detector::Union{Detector, Vector{Detector}},
    mc::AbstractArray,
    eta::AbstractArray,
    chi1::AbstractArray,
    chi2::AbstractArray,
    dL::AbstractArray,
    theta::AbstractArray,
    phi::AbstractArray,
    iota::AbstractArray,
    psi::AbstractArray,
    tcoal::AbstractArray,
    phiCoal::AbstractArray,
    Lambda1=nothing,
    Lambda2=nothing;
    fmin=2.0,
    fmax = nothing,
    res = 1000,
    useEarthMotion = false,
    rho_thres=12.,
    alpha = 0.0,
    coordinate_shift = true,
    return_SNR = false,
    auto_save =false,
    name_folder = nothing,
    save_catalog = false,
)
    nEvents = length(mc)    

    if name_folder === nothing
        name_folder = _event_type(model) 
    end

    if Lambda1===nothing
        Lambda1 = zeros(size(mc))
    end
    if Lambda2===nothing
        Lambda2 = zeros(size(mc))
    end

    if model == TaylorF2()
        nPar = _npar(model, Lambda1[1], Lambda2[1])
    else
        nPar = _npar(model)
    end

    Fishers = Array{Float64}(undef, nEvents, nPar, nPar)
    if return_SNR == true
        SNRs = Array{Float64}(undef, nEvents)
        elapsed_time = @elapsed @showprogress desc="Computing Fishers and SNRs..." @threads for ii in 1:nEvents  
            Fishers[ii,:,:], SNRs[ii] = FisherMatrix(
                model,
                detector, 
                mc[ii], 
                eta[ii], 
                chi1[ii], 
                chi2[ii], 
                dL[ii], 
                theta[ii], 
                phi[ii], 
                iota[ii], 
                psi[ii], 
                tcoal[ii], 
                phiCoal[ii], 
                Lambda1[ii], 
                Lambda2[ii], 
                fmin=fmin, 
                fmax=fmax, 
                res = res, 
                useEarthMotion = useEarthMotion, 
                rho_thres=rho_thres, 
                alpha = alpha, 
                coordinate_shift = coordinate_shift, 
                return_SNR=true
            )
        end 

        println("Fisher matrices and SNRs computed!")
        if elapsed_time > 60.0
            elapsed_time = elapsed_time/60
            println("The evaluation took: ", elapsed_time, " minutes.")
        else
            println("The evaluation took: ", elapsed_time, " seconds.")
        end

        if auto_save == true  
            path = pwd()
            mkpath("output/"*name_folder)
            path = pwd()*"/output/"*name_folder*"/"
            date = Dates.now()
            date_format = string(Dates.format(date, "e dd u yyyy HH:MM:SS"))
            h5open(path*"Fishers_SNRs.h5", "w") do file
                write(file, "Fishers", Fishers)
                attributes(file)["number_events"] = nEvents
                if typeof(detector) == Vector{Detector}
                    label = [detector[i].label for i in eachindex(detector)]
                    attributes(file)["Detectors"] = label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers and the SNRs for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*join(label, ",")*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format

                else
                    attributes(file)["Detectors"] = detector.label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers and the SNRs for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*detector.label*" detector and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format

                end     
                write(file, "SNRs", SNRs) 
                if save_catalog 
                    write(file, "mc", mc)
                    write(file, "eta", eta)
                    write(file, "chi1", chi1)
                    write(file, "chi2", chi2)
                    write(file, "dL", dL)
                    write(file, "theta", theta)
                    write(file, "phi", phi)
                    write(file, "iota", iota)
                    write(file, "psi", psi)
                    write(file, "tcoal", tcoal)
                    write(file, "Lambda1", Lambda1)
                    write(file, "Lambda2", Lambda2)
                end
            end
        end
        return Fishers, SNRs
    else
        elapsed_time = @elapsed  @showprogress desc="Computing Fishers..."  @threads for ii in 1:nEvents  
                    Fishers[ii,:,:]=FisherMatrix(
                        model,
                        detector,
                        mc[ii],
                        eta[ii], 
                        chi1[ii],
                        chi2[ii],
                        dL[ii],
                        theta[ii],
                        phi[ii], 
                        iota[ii], 
                        psi[ii], 
                        tcoal[ii], 
                        phiCoal[ii],
                        Lambda1[ii],
                        Lambda2[ii],
                        fmin=fmin, 
                        fmax=fmax, 
                        res = res, 
                        useEarthMotion = useEarthMotion, 
                        rho_thres=rho_thres, 
                        alpha = alpha, 
                        coordinate_shift = coordinate_shift
                    )
                end 
        println("Fisher matrices computed!")
        if elapsed_time > 60.0
            elapsed_time = elapsed_time/60
            println("The evaluation took: ", elapsed_time, " minutes.")
        else
            println("The evaluation took: ", elapsed_time, " seconds.")
        end
        if auto_save == true  
            path = pwd()
            mkpath("output/"*name_folder)
            path = pwd()*"/output/"*name_folder*"/"
            date = Dates.now()
            date_format = string(Dates.format(date, "e dd u yyyy HH:MM:SS"))
            h5open(path*"Fishers.h5", "w") do file
                write(file, "Fishers", Fishers)
                attributes(file)["number_events"] = nEvents
                if typeof(detector) == Vector{Detector}
                    label = [detector[i].label for i in eachindex(detector)]
                    attributes(file)["Detectors"] = label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*join(label,",")*" detectors and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format
                else
                    attributes(file)["Detectors"] = detector.label
                    attributes(file)["What_this_file_contains"] = "This file contains the Fishers for "*string(nEvents)*" events obtained with the "*string(typeof(model))*" waveform model. The calculations are performed with the "*detector.label*" detector and the correction due to Earth Motion was "*string(useEarthMotion)*"."
                    attributes(file)["date"] = date_format
                end    
                if save_catalog
                    write(file, "mc", mc)
                    write(file, "eta", eta)
                    write(file, "chi1", chi1)
                    write(file, "chi2", chi2)
                    write(file, "dL", dL)
                    write(file, "theta", theta)
                    write(file, "phi", phi)
                    write(file, "iota", iota)
                    write(file, "psi", psi)
                    write(file, "tcoal", tcoal)
                    write(file, "Lambda1", Lambda1)
                    write(file, "Lambda2", Lambda2)
                end
            end
        end
        return Fishers
    end
        
end


"""
This function reads the Fisher matrices and/or the SNRs
    _read_Fishers_SNRs(path, SNR=true)

    #### Input arguments:
    -  `path` : string, path to the file
    -  `SNR` : bool, default true, if true the function reads also the SNRs

    #### Output:
    - `Fishers`  : matrix, Fisher Matrix
    - `SNRs`  : vector, SNRs (if SNR=true)

    #### Example:
    ```julia
    Fishers, SNRs = _read_Fishers_SNRs("output/BBH/Fishers_SNRs.h5")
    ```

"""
function _read_Fishers_SNRs(path; SNR=true)
    if SNR == true
        Fishers, SNRs = h5open(path, "r") do file
            println("Attributes: ", keys(attributes(file)))
            attributes_keys = keys(attributes(file))
            for i in 1:length(keys(attributes(file)))
                println(attributes_keys[i], ": ", read(attributes(file)[attributes_keys[i]]))
            end
            println("Keys: ", keys(file))
            Fishers = read(file, "Fishers")
            SNRs = read(file, "SNRs")
            return Fishers, SNRs
        end
    else
        Fishers = h5open(path, "r") do file
            println("Attributes: ", keys(attributes(file)))
            attributes_keys = keys(attributes(file))
            for i in 1:length(keys(attributes(file)))
                println(attributes_keys[i], ": ", read(attributes(file)[attributes_keys[i]]))
            end
            println("Keys: ", keys(file))
            Fishers = read(file, "Fishers")
            return Fishers
        end
    end
end





end# of module
