from __future__ import division
import numpy as np


def MTF(spatial_frequency,eccentricity, paper='Williams1996_astig'):
    """    
    Compute the modulation frequency transfer as a function of eccentricity
    based on the equation derrived in Navarro, Artal, and Williams 1993. 
    
    .. math::
       F = 1-C * exp(-A*f) + C*exp(-B*f)
    
    :param spatial_frequency: array or integer of spatial frequencies
    :param eccentricity: eccentricity from 0 to 60 degrees at which to compute MTF.\
    Linear interpolation is used for values that fall inbetween those \
    reported in Navarro et al and are, therefore, less reliable though as expected. \
    Eccentricities reported (in degrees): [0, 5, 10, 20, 30, 40, 50, 60]
    :param paper: Choose the coefficients to use for generating the curve.  \n
                 Options are : \n
                 * 'Williams1996_astig' \n
                 * 'Williams1996_clc' \n
                 * 'Navarro1993' 
                 
    :returns: MTF for input frequencies at given eccentricity.

    
    """
    if paper == 'Navarro1993':    
        theta = [0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0] # in degrees
        A_theta = [0.172, 0.245, 0.245, 0.328, 0.606, 0.82, 0.93, 1.89] # in degrees
        B_theta = [0.037, 0.041, 0.041, 0.038, 0.064, 0.064, 0.059, 0.108] # in degrees
        C_theta = [0.22, 0.2, 0.2, 0.14, 0.12, 0.09, 0.067, 0.05] # unitless

    if paper == 'Williams1996_astig':
        theta = [0.0, 10.0, 20.0,  40.0]
        A_theta = [0.0129, 0.0140, 0.0082, 0.0059]
        B_theta = [0.0816, 0.1036, 0.1313, 0.1555]
        C_theta = [0.7921, 0.8299, 0.9120, 0.9178]
 
    if paper == 'Williams1996_clc':
        theta =   [0.0,    10.0,   20.0,   40.0]
        A_theta = [0.0122, 0.0154, 0.0000, 0.0000]
        B_theta = [0.0988, 0.1466, 0.2305, 0.4663]
        C_theta = [0.8172, 0.8266, 0.9378, 0.9515]       
    # linear interpolation:

    A = np.interp(eccentricity,theta,A_theta)
    B = np.interp(eccentricity,theta,B_theta)
    C = np.interp(eccentricity,theta,C_theta)

    f = (1.0 - C) * np.exp(-A * spatial_frequency) + C * np.exp(-B * spatial_frequency)
    
    return f
    

    
def MTF_Pupil(spatial_frequency, pupilDiameter, normalized = 0):
    """
    
    Compute the modulation frequency transfer as a function of pupil size
    as derrived in Artal and Navarro 1994. 
    F = 1-C * exp(-A*f) + C*exp(-B*f)

    :param spatial_frequency: array or integer of spatial frequencies
    :param eccentricity: eccentricity from 0 to 60 degrees at which to compute MTF.\
    Linear interpolation is used for values that fall inbetween those \
    reported in Navarro et al and are, therefore, less reliable though as expected. \
    Pupil size reported (in mm): [2.5, 3.0, 4.0, 6.0, 8.0]. The behavior of \
    this function in between reported sizes may be worse than MTF(). \
    :param normalized: Default = 0 (no). 1 (yes) will normalize the spatial \
    frequency.  
    
    .. note::
       * This is not fully tested.  Permits comparison of a given
         optical system (eye) to a perfect one.  Requires consideration of maximum 
         spatial frequency of the optical system (uLim).    
    

    :returns f: MTF for input frequencies at given eccentricity.
    
    .. warning:: 
       normalized option 1 is not fully tested.
    
    """

    pupil = [2.5, 3.0, 4.0, 6.0, 8.0] # in mm.

    if normalized == 0:

        A_pupil = [0.16, 0.16, 0.18, 0.31, 0.53] # in degrees
        B_pupil = [0.06, 0.05, 0.04, 0.06, 0.08] # in degrees
        C_pupil = [0.36, 0.28, 0.18, 0.2, 0.11] # unitless
        
        A = np.interp(pupilDiameter,pupil,A_pupil)
        B = np.interp(pupilDiameter,pupil,B_pupil)
        C = np.interp(pupilDiameter,pupil,C_pupil)

        f = (1 - C) * np.exp(-A * spatial_frequency) + C * np.exp(-B * spatial_frequency)

    elif normalized == 1:
        
        A_pupil = [10.57, 12.68, 19.04, 49.19, 112.15] # no dimension
        B_pupil = [3.96, 3.96, 4.23, 9.52, 16.92] # no dimension
        C_pupil = [0.36, 0.28, 0.18, 0.2, 0.11] # no dimension
        uLim = [66.1, 79.3, 105.8, 158.7, 211.6] # in cycles/degree;
        
        A = np.interp(pupilDiameter,pupil,A_pupil)
        B = np.interp(pupilDiameter,pupil,B_pupil)
        C = np.interp(pupilDiameter,pupil,C_pupil)
        u = np.interp(pupilDiameter,pupil,uLim)
    
        f = (1 - C) * np.exp(-A * (spatial_frequency / u)) + C * np.exp(-B * (spatial_frequency / u))
    
    return f



def spectsens(LambdaMax = 559, OpticalDensity = 0.2, Output = 'log', StartWavelength = 380, 
                EndWavelength = 780, Res = 1000):
    """
    
    This function returns a photopigment spectral sensitivity curve as defined by 
    Carroll, McMahon, Neitz, and Neitz. 

    
    :param LambdaMax: Wavelength peak for photopigment
    :param OpticalDensity: optical density required
    :param OutputType:  log or anti-log.  \n
                        if log, maximum data ouput is 0. \n
                        if anti-log, data output is between 0 and 1. \n                   
    :param StartWavelength: beginning wavelength
    :param EndWavelength: end wavelength
    :param Resolution: Number of data points
    
    :returns: array of sensitivity values.
    :rtype: np.array
    
    .. note::
       Ported from Jim K's Matlab function.
       
    """

    A = 0.417050601
    B = 0.002072146
    C = 0.000163888
    D = -1.922880605
    E = -16.05774461
    F = 0.001575426
    G = 5.11376E-05
    H = 0.00157981
    I = 6.58428E-05
    J = 6.68402E-05
    K = 0.002310442
    L = 7.31313E-05
    M = 1.86269E-05
    N = 0.002008124
    O = 5.40717E-05
    P = 5.14736E-06
    Q = 0.001455413
    R = 4.217640000E-05
    S = 4.800000000E-06
    T = 0.001809022
    U = 3.86677000E-05
    V = 2.99000000E-05
    W = 0.001757315
    X = 1.47344000E-05
    Y = 1.51000000E-05
    Z = OpticalDensity+0.00000001

    
    A2 = (np.log10(1.0 / LambdaMax) - np.log10(1.0 / 558.5))
    
    vector = np.log10(np.linspace(StartWavelength,EndWavelength,Res)**-1.0)
    
    const = 1.0 / np.sqrt(2.0 * np.pi)

    exTemp = (np.log10(-E + E * np.tanh(-(((10.0**(vector - A2))) - F) / G)) + D + 
                        A * np.tanh(-(((10.0**(vector - A2))) - B) / C) - 
                        (J / I * (const * np.exp(1.0)**(-0.5 * (((10.0**(vector - A2)) - H) / I)**2.0))) - 
                        (M / L * (const * np.exp(1.0)**(-0.5 * (((10.0**(vector - A2)) - K) / L)**2.0))) - 
                        (P / O * (const * np.exp(1.0)**(-0.5 * (((10.0**(vector - A2)) - N) / O)**2.0))) + 
                        (S / R * (const * np.exp(1.0)**(-0.5 * (((10.0**(vector - A2)) - Q) / R)**2.0))) + 
                        ((V / U * (const * np.exp(1.0)**(-0.5 * (((10.0**(vector - A2)) - T) / U)**2.0))) / 10.0) + 
                        ((Y / X * (const * np.exp(1.0)**(-0.5 * (((10.0**(vector - A2)) - W) / X)**2.0))) / 100.0))
    ODTemp = np.log10((1.0 - 10.0** -((10.0**exTemp) * Z)) / (1.0 - 10**-Z))

    if Output.lower() == 'log':
        extinction = exTemp
        withOD = ODTemp
    else:
        extinction = 10.0**(exTemp)
        withOD = 10.**(ODTemp)
    
    return withOD[0], extinction
    

def StilesCrawford1stKind(x, xmax, n, rho):
    """
    
    Stiles-Crawford effect of the 1st order.
    
    :param x: array of x values
    :param xmax: maximum value
    :param n: 
    :param rho: angle
    
    :returns: styles-crawford effect
    
    .. warning::
       This funciton is not finished or tested.
    """
    
    return np.log(n) - rho*(x - xmax)**2

    
    
def SterhlRatio(sigma, lam):
    """
    Find the Sterhl ratio.
    
    .. warning::
       This funciton is not finished or tested.
       
    """
    
    return np.exp(-2.0 * np.pi * sigma / lam)**2.0

    

def NumericalAperature(n, theta):
    """
    Find the numerical aperature of a system
    
    .. warning::
       Not really working. Untested.
       
    """
    return n * np.sin(theta)

    
def DiffractionLimit(lam, n,theta):
    """
    
    Find the diffraction limit of a system.
    
    :param lam: wavelength
    :param n: 
    :param theta:
    
    .. warning::
       This funciton is not finished or tested.
        
    """
    return lam / 2.0 * ( NumericalAperature(n, theta) )
    

