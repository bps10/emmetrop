from __future__ import division
import numpy as np

from base import optics as o
from base.data import rad2deg
from base import cones as cones

# emmetrop imports:
from emmetrop.analysis import Information as info
from emmetrop.eye.eyeModel import traceEye
from emmetrop.eye.movement import brownian_motion
from emmetrop.eye.filter import gauss

def genMeta():
    '''Just temporary.
    '''
    # get from ray tracer eventually
    _meta = {}
    _meta['samples'] = 399
    _meta['retImg'] = 0.1995 #this is only radius, D = 2x 0.1995
    _meta['eye_length'] = 23.92
    _meta['pupil_size'] = 3
    radians = 2 * np.arctan(_meta['retImg'] / _meta['eye_length'])
    _meta['deg'] = rad2deg(radians)
    _meta['mm/deg'] = 2 * _meta['retImg'] / _meta['deg']

    cycles = np.arange(0, _meta['samples']) / 2.0
    cpd = (cycles / _meta['retImg']) * _meta['mm/deg'] 

    return _meta, cpd

def NeitzModel(ImageData, rec_field, cpd, field_angle, _meta, analysis_args, glasses):
    '''This function organizes the entire operation.
    A dictionary self.Analysis is created to reflect the 
    user options. All subsequent methods will use the keys
    of this dictionary.
    '''
    Analysis = {}

    if 'dist' in analysis_args:
        dist_range = 10 ** (np.arange(5, 23, 0.5) / 3.0)
    else:
        dist_range = np.array([1e8])

    if 'focus' in analysis_args:
        focus_range = np.arange(0, 3.6, 0.5)
    else: 
        focus_range = np.array([0])

    if 'off_axis' in analysis_args:
        axis_range = np.arange(0, 21, 2)
    else:
        axis_range = np.array([field_angle])

    if 'pupil_size' in analysis_args:
        pupil_range = np.arange(2, 9, 1)
    else:
        pupil_range = np.array([3])

    wavelen = rec_field.keys()
    wavelen.remove('max')
    wavelen = np.sort(wavelen)

    j = 0
    for dist in dist_range:
        for focus in focus_range:
            for axis in axis_range:
                for pupil in pupil_range:
                    for wave in wavelen:

                        Analysis[j] = {
                            'dist': dist,
                            'focus': focus,
                            'off_axis': axis,
                            'pupil_size': pupil, 
                            'wavelength': wave,
                            'line': addLineStyle(dist, focus, axis, pupil), }
                        j += 1

                        # generate new rec field base on wavelen if necessary

    # then find biggest response, normalize everything else to that
    
    Analysis, diffract = computeConeActivity(Analysis, ImageData,
        rec_field, cpd, _meta, True, glasses)      
    Analysis = totalActivity(Analysis, diffract)
    Analysis = estimateInfo(Analysis, ImageData, diffract)

    return Analysis, diffract

def addLineStyle(dist, focus, axis, pupil):
    """Add line style for use with plots
    """
    r = 0 #focus / 2
    g = 0 #np.log10(dist) / (25 / 3)
    b = 0 #axis / 20
    a = 0.4
    rgb = [r, g, b, a]
    line = {'style': '-', 'color': rgb}
    return line
     

def computeConeActivity(Analysis, ImageData, rec_field, cpd, _meta,
                _brownian=True, glasses=False):
    """Compute the estimated activity of a cone photoreceptor.
    
    :param Receptive_Field: a handle to the spline fitted receptive field.
    :type Receptive_Field: function handle.
    """
    Rec_Field = (rec_field[540]['fft'] / rec_field['max']) 
        
    ImageData['fitLaw'] = ImageData['powerlaw'](cpd[1:])
    powerlaw = ImageData['fitLaw']

    if _brownian:    
        temp = np.arange(1, 80)
        movement_filter = brownian_motion(cpd[1:], temp)
        powerlaw *= movement_filter

    # compute the diffraction limited case seperately
    diffract = {}
    diff, _x = o.diffraction(_meta['samples'], 
                                    _meta['pupil_size'],
                                    16.6, 
                                    ref_index=1.4, 
                                    wavelength=550.0)
    # now interpolate mtf into cpd's of analysis:
    mtf = np.interp(cpd, _x, diff)

    # remove zeros to avoid errors in future computations:
    ind = np.where(mtf != 0)[0]
    diffract['cpd'] = cpd[ind]
    diffract['mtf'] = mtf[ind]

    diffract['preCone'] =  (powerlaw[ind] * 
            diffract['mtf'][ind])
    diffract['retina'] = (diffract['preCone'] *
                                    Rec_Field[ind])

    if glasses:
        # if glasses on, compute effect after diffraction case
        powerlaw *= gauss(cpd[1:], 10)

    for key in Analysis:
        # find cone fft:
        wv = Analysis[key]['wavelength']
        Rec_Field = (rec_field[wv]['fft'] / rec_field['max'])

        # generate MTFs for each condition:
        intensity = traceEye(
                            Analysis[key]['dist'], 
                            Analysis[key]['off_axis'], 
                            Analysis[key]['pupil_size'], 
                            Analysis[key]['focus'],
                            Analysis[key]['wavelength'])
        psf = o.genPSF(intensity, _meta['samples'])[1]
        Analysis[key]['mtf'] = o.genMTF(psf)[ind]

        Analysis[key]['preCone'] = (powerlaw[ind] * 
            Analysis[key]['mtf'])

        Analysis[key]['retina'] = (Analysis[key]['preCone'] *
                                    Rec_Field[ind])

    return Analysis, diffract
                  
def totalActivity(Analysis, diffract, print_opt=True):
    """Compute the estimated activity in a photoreceptor.
    
    :param self: no params yet.
    
    :returns: percent of integrated.
    
    .. note:: 
       Completely hard coded right now.  Would like to change that.
    
    """
    diffraction_total = np.sum(diffract['retina'])
    print diffraction_total
    for key in Analysis:
        total = np.sum(Analysis[key]['retina'])
        Analysis[key]['percent'] = (total /
                                diffraction_total)

    if print_opt:
        print ' '
        print 'Photoreceptor Activity (proportion of diffraction)'
        print ' '
        print 'index\tobj_dist_mm\tfocus_D\toff_axis_deg\twavelen_nm\tpupil_mm\tproportion'
        for key in Analysis:
            line = (str(key) + 
                '\t' + str(round(Analysis[key]['dist'], 3)) +
                '\t' + str(round(Analysis[key]['focus'], 3)) +
                '\t' + str(round(Analysis[key]['off_axis'], 3)) +
                '\t' + str(round(Analysis[key]['wavelength'], 3)) +
                '\t' + str(round(Analysis[key]['pupil_size'],1)) +
                '\t' + str(round(Analysis[key]['percent'], 4)))
            print line

    return Analysis

        
def estimateInfo(Analysis, ImageData, diffract, print_opt=False):
    """Estimate the information in a simple linear cone receptive field.
    
    Information is estimated with Garrigan et al.'s Gaussian approximation
    method.
    
    .. math::
       I_1(S,E) = \\frac{1}{2}\log_2(1+SNR) , 
    
    where :math:`SNR` is the ratio of signal power to noise power.
    
    Information in an array of cones is then represented:
    
    .. math::
       I_N(S,E) = I_1(S,E)N^\delta ,
    
    with :math:`\\delta` representing a derived scale factor.
    
    :param Receptive_Field: type of receptive field to use (FFT, Jay)
    :param print_opt: decide whether to print results (True) or not (False)
    :type print_opt: bool

    
    :returns: plot of information estimate
    
    .. warning::
       This function is under development. 
    
    :usage: estimateInfo is called by ComputeConeActivity() function.

    **This function produces:**
    
    .. figure:: ../../Figures/InfoPlot.png 
       :height: 300px
       :width: 400px
       :align: center   
       
    """        
    total_images = ImageData['totalImages']        

    for key in Analysis:
        Analysis[key]['cones'] = [1,2,3,4,5,6]
        Analysis[key]['info'] = np.zeros(len(
                                            Analysis[key]['cones']))

    for amp in ImageData['rawAmp']:

        for key in Analysis:
            ind = len(diffract['cpd'])
            fooInfo = info.SingleConeEntropyFunc((amp[ind] ** 2 *
                                        Analysis[key]['retina']), 
                                        Analysis[key]['cones'])      
            Analysis[key]['info'] += fooInfo / total_images

    if print_opt == True:
        print ' '
        print 'Information'
        print '------------'
        for key in Analysis:
            print key, ': ', Analysis[key]['info']

    return Analysis


if __name__ == "__main__":

    eye = SchematicEyeAnalysis()
