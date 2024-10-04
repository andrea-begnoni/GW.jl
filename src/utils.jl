module UtilsAndConstants    # this is the name of the module

using ForwardDiff # these are the Julia packages that are used in this module
using LinearAlgebra
using LaTeXStrings


export GMsun_over_c3, GMsun_over_c2, uGpc, GMsun_over_c2_Gpc, REarth_km, clight_kms, clightGpc, Lamt_delLam_from_Lam12,
         _ra_dec_from_theta_phi_rad, _theta_phi_from_ra_dec_rad, CovMatrix, Errors, SkyArea, _orientationBigCircle


##############################################################################
GMsun_over_c3 = 4.925491025543575903411922162094833998e-6 # seconds

GMsun_over_c2 = 1.476625061404649406193430731479084713e3 # meters

uGpc = 3.085677581491367278913937957796471611e25 # meters

GMsun_over_c2_Gpc = GMsun_over_c2 / uGpc # Gpc

REarth_km = 6371.00 # km

clight_kms = 2.99792458e5 # km/s

clightGpc = clight_kms / 3.0856778570831e+22

"""Solar mass"""
MSUN = 1.988409902147041637325262574352366540e30  # kg

"""Geometrized nominal solar mass, m"""
MRSUN = GMsun_over_c2


##############################################################################


r"""
Compute the dimensionless tidal deformability combinations L"\tilde{Lambda}" and L"\delta\tilde{Lambda}", defined in `arXiv:1402.5156 <https://arxiv.org/abs/1402.5156>`_ eq. (5) and (6), as a function of the dimensionless tidal deformabilities of the two objects and the symmetric mass ratio.
Function from GWFAST.

#### Input argument:
-  float Lambda1 Tidal deformability of object 1
-  float Lambda2 Tidal deformability of object 2
-  float eta The symmetric mass ratio of the object.

#### Outputs:
- (float, float) L"\tilde{Lambda}" and L"\delta\tilde{Lambda}", 

"""
function Lamt_delLam_from_Lam12(Lambda1, Lambda2, eta)

    eta2 = eta*eta

    Seta = sqrt(1.0 - 4.0 * eta)  
        
    Lamt = (8. /13.)*((1. + 7. *eta - 31. *eta2)*(Lambda1 + Lambda2) + Seta*(1. + 9. *eta - 11. *eta2)*(Lambda1 - Lambda2))
    
    delLam = 0.5*(Seta*(1. - 13272. /1319. *eta + 8944. /1319. *eta2)*(Lambda1 + Lambda2) + (1. - 15910. /1319. *eta + 32850. /1319. *eta2 + 3380. /1319. *eta2*eta)*(Lambda1 - Lambda2))
    
    return Lamt, delLam

end




"""
Function to obtain RA and DEC from theta and phi in radians. Result is in radians.

#### Input arguments:
-  float theta Polar angle in radians.
-  float phi Azimuthal angle in radians.

#### Outputs:
-  (float, float) RA and DEC in radians.

#### Example:
```julia
ra, dec = _ra_dec_from_theta_phi_rad(0.1, 0.2)
```

"""
function _ra_dec_from_theta_phi_rad(theta::Union{Float64, ForwardDiff.Dual}, phi::Union{Float64, ForwardDiff.Dual})

    ra = phi
    dec = 0.5 * pi - theta
    return ra, dec
end

function _ra_dec_from_theta_phi_rad(theta::AbstractArray, phi::AbstractArray)
    ra = phi
    dec = [0.5 * pi - theta[i] for i in eachindex(theta)]
    return ra, dec
end

function _theta_phi_from_ra_dec_rad(ra::Union{Float64, ForwardDiff.Dual}, dec::Union{Float64, ForwardDiff.Dual})
    theta = 0.5 * pi - dec
    phi = ra
    return theta, phi
end

function _theta_phi_from_ra_dec_rad(ra::AbstractArray, dec::AbstractArray)
    theta = [0.5 * pi - dec[i] for i in eachindex(dec)]
    phi = ra
    return theta, phi
end

"""
Function to compute the covariance matrix from the Fisher matrix. If the matrix is not positive definite, it returns a zero matrix. This can happen when the Fisher matrix is not invertible (e.g., when the signal-to-noise ratio is low).
If the Fisher matrix is made of zeros, it returns a zero matrix.
#### Input arguments:
-  Fisher matrix: square matrix of size NxN.

#### Outputs:
-  Covariance matrix: square matrix of size NxN.

#### Example:
```julia
    covMatrix = CovMatrix([1 0; 0 1])
```

"""
function CovMatrix(Fisher::Matrix{Float64}, print_error = true)
    covMatrix = zeros(size(Fisher))  # Initialize covMatrix as a zero matrix of the same size as Fisher

    if Fisher[1,1]==0.
        return covMatrix
    end
    # normalize the fisher, technique by GWFISH 
    diagonal = diag(Fisher)
    Fisher = Fisher ./ sqrt.(diagonal * diagonal')
    try
        covMatrix = inv(cholesky(Fisher))
    catch
        if print_error == true
            println("Fisher matrix is not positive definite. Returning zero matrix.")
        end
    end
    covMatrix = covMatrix ./ sqrt.(diagonal * diagonal')

    return covMatrix
end

"""
Same as CovMatrix(Fisher::Matrix{Float64}) but for a 3D array of Fisher matrices.
"""
function CovMatrix(Fisher::Array{Float64, 3}, print_error = false)
    covMatrix = zeros(size(Fisher))  # Initialize covMatrix as a zero matrix of the same size as Fisher

    n = size(Fisher)[1]
    for i in 1:n
        covMatrix[i,:,:] = CovMatrix(Fisher[i,:,:], print_error)
    end
    
    return covMatrix
end
"""
Calculate the errors from the covariance matrix. The errors are the square root of the diagonal elements of the covariance matrix.

#### Input arguments:
-  covMatrix: square matrix of size NxN.

#### Outputs:
-  errors: array of size N.

#### Example:
```julia
    errors = Errors([1 0; 0 4]) # returns [1.0, 2.0]
"""
function Errors(covMatrix::Matrix{Float64})
    return sqrt.(diag(covMatrix))
end

"""
Same as Errors(covMatrix::Matrix{Float64}) but for a 3D array of covariance matrices.
"""
function Errors(covMatrix::Array{Float64, 3})
    
    nCov = size(covMatrix)[1]
    nPar = size(covMatrix)[2]
    errors = zeros(nCov, nPar)

    for ii in 1:nCov
        errors[ii,:] = sqrt.(diag(covMatrix[ii,:,:]))
    end

    return errors
end

"""
Calculate the sky area from the covariance matrix and the catalog of sources. The sky area is the area of the error ellipse in the sky, in square degrees.

#### Input arguments:
-  covMatrix: square matrix of size NxN.
-  thetaCatalog: float with the polar angle of the source in the catalog.

#### Optional arguments:
-  percent_level: float with the confidence level. Default is 90%.

#### Outputs:
-  float with the sky area in square degrees.

#### Example:
```julia
    covMatrix = Matrix{Float64}(I, 11,11)
    sky_area = SkyArea(covMatrix, 0.1) 
"""
function SkyArea(covMatrix::Matrix{Float64}, thetaCatalog; percent_level = 90)
    # take phi and theta from covMatrix and multiply them

    sigmaTheta_sq  = covMatrix[6,6]
    sigmaPhi_sq = covMatrix[7,7]
    sigmaThetaPhi = covMatrix[6,7]
    
    # From Barak, Cutler, PRD 69, 082005 (2004), gr-qc/0310125
    # results in grad
    return - 2*pi*sqrt(sigmaTheta_sq*sigmaPhi_sq - sigmaThetaPhi^2)*abs(sin(thetaCatalog))*log(1-percent_level/100)*(180/pi)^2
end

function SkyArea(covMatrix::Array{Float64,3}, thetaCatalog::Array; percent_level = 90)
    
    nCov = size(covMatrix)[1]
    skyArea = zeros(nCov)

    for ii in 1:nCov
        skyArea[ii] = SkyArea(covMatrix[ii,:,:], thetaCatalog[ii], percent_level=percent_level)
    end

    return skyArea
end

""" 
Calculate the optimal angle for the orientation of the detector to maximize the signal from CBC gravitational wave. 
The function can be used in the definition of 2 detectors network.
It outputs the orientation of the second detectors such that it is at 45 degrees 
These reasoning works if the first detector is oriented towards the East.
If you want to use it for a different orientation, you need to add the angle of the first detector to the output of this function.

#### Input arguments:
- DetectorStructure det1: DetectorStructure object of the first detector.
- DetectorStructure det2: DetectorStructure object of the second detector.

#### Outputs:
-  float orientation_CBC Optimal angle for the optimal orientation of the second detector in radians

#### Example:
```julia
    orientation_CBC = _orientationBigCircle(0.1, 0.2, 0.3, 0.4)
 ```
""" 
function _orientationBigCircle(det1, det2)

    lat1 = det1.latitude_rad
    long1 = det1.longitude_rad
    orientation1 = det1.orientation_rad
    lat2 = det2.latitude_rad
    long2 = det2.longitude_rad


    equatorialRadiusEarth = 6378.137 #from https://nssdc.gsfc.nasa.gov/planetary/factsheet/earthfact.html 
    polarRadiusEarth = 6356.752 

    eFactor = 1 - (equatorialRadiusEarth^2 - polarRadiusEarth^2) / equatorialRadiusEarth^2  # 1-e^2
    eS = (equatorialRadiusEarth^2 - polarRadiusEarth^2) / equatorialRadiusEarth^2 # e^2

    Lambda1 = eFactor * tan(lat2)/tan(lat1) + eS*sqrt( (1 + eFactor * tan(lat2)^2) / (1 + eFactor * tan(lat1)^2))
    Lambda2 = eFactor * tan(lat1)/tan(lat2) + eS*sqrt( (1 + eFactor * tan(lat1)^2) / (1 + eFactor * tan(lat2)^2))


    # alpha1 = atan( sin(long2 - long1) , ((Lambda1 - cos(long2 - long1))*sin(lat1)))
    # alpha2 = atan( sin(long1 - long2) , ((Lambda2 - cos(long1 - long2))*sin(lat2)))
    alpha1 = atan( sin(long2 - long1) , ((Lambda1 - cos(long2 - long1))*sin(lat1)))
    alpha2 = atan( sin(long1 - long2) , ((Lambda2 - cos(long1 - long2))*sin(lat2)))

    return mod(orientation1 + alpha1 - alpha2 + pi/4, pi)

end


end