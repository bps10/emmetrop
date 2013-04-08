import numpy as np
from scipy import interpolate


class ConeReceptiveFields(object):
    """This class generates models of cone receptive fields.
    
    Currently there is no real need for a class. However, it is expected
    that as this project grows, we will introduce additional models and this
    class will surve a more dynamic purpose.
    
    if no cone_spacing or locations passed, defaults::
       cone_spacing=[2.0, 0.5]
       location = ['periph', 'fovea']
    """
    def __init__(self, freqs, cone_spacing=[2.0, 0.5],
                 location = None):
        """
        """
        if len(freqs) != len(cone_spacing):
            raise ValueError("Sorry, length of location must match \
            cone_spacing")
        if not location:
            location = ['periph', 'fovea']
            
        self.genReceptiveFields(freqs, cone_spacing, location)
            
    def genReceptiveFields(self, freqs, cone_spacing, location):
        """Create a difference of gaussians receptive field and plot it.
    
        :param excite_SD: standard deviation parameter of excitatory gaussian.
        :param inhibit_SD: standard deviation parameter of inhibitory \
        surround gaussian. Default = 5.0.   
        
        :returns: RF_DOG, RF_SPLINE, FFT_RF

        """
        
        # traditional RF
        N = 300
        Xvals = np.arange(-15, 15, 10. / (2.0 * N))
    
        RF_DOG = {}
        FFT_RF = {}
        RF_SPLINE = {}
        RField = {}
    
        RF_DOG = findSpacing(Xvals, cone_spacing)
        FFT_RF = Fourier(RF_DOG, N)
        length = np.floor(FFT_RF.shape[0] / 2.) + 1
        RField = normRField(freqs, Xvals[length:] * 60,
                                 FFT_RF[length:])

        for i, loc in enumerate(location):
            
            RF_DOG[loc] = findSpacing(Xvals, cone_spacing[i])
            FFT_RF[loc] = Fourier(RF_DOG[loc], N)
            length = np.floor(FFT_RF[loc].shape[0] / 2.) + 1 
            RField[loc] = normRField(freqs, Xvals[length:] * 60,
                                     FFT_RF[loc][length:])

            ## spline interpolation handle            
            RF_SPLINE[loc] = interpolate.splrep(Xvals[length:] * 60,
                                            FFT_RF[loc][length:], s=0)

        # manual method --- Jay's method:        
        spatial_frequencies = np.arange(1, length - 1)
    
        # allocate memory        
        sine_wave = np.zeros((Xvals.shape[0], spatial_frequencies.shape[0]))
        Jay_RF = {}
        Jay_RField = {}
        Jay_CR = {}
        for i, loc in enumerate(location):
            cone_response = np.zeros(spatial_frequencies.shape[0])        
            for i, thisFrequency in enumerate(spatial_frequencies):
                    
                sine_wave[:,i] = sineWave(thisFrequency, Xvals)
                
                cone_response[i] = np.sum(sine_wave[:,i] * RF_DOG[loc])

            Jay_CR[loc] = cone_response / np.sum(cone_response)
            Jay_RF[loc] = interpolate.splrep(Xvals[length:] * 60,
                                                    cone_response, s=0)

            Jay_RField[loc] = normJayRField(freqs, Jay_RF[loc])

            
    
        self.receptive_field =  {
                            'length': length,
                            'dog': RF_DOG,
                            'coneResponse': {'fft': FFT_RF, 'jay': Jay_CR},
                            'spline': {'fft':RF_SPLINE, 'jay': Jay_RF},
                            'RField': {'fft': RField, 'jay': Jay_RField},
                            'sineWave': sine_wave,
                            'xvals': Xvals
                            }
                            
    def returnReceptiveField(self):  
        """
        Return a dictionary of receptive field data.
        """                  
        return self.receptive_field



def genReceptiveFields(freqs, cone_spacing):
    """Create a difference of gaussians receptive field and plot it.

    :param excite_SD: standard deviation parameter of excitatory gaussian.
    :param inhibit_SD: standard deviation parameter of inhibitory \
    surround gaussian. Default = 5.0.

    :returns: RF_DOG, RF_SPLINE, FFT_RF

    """

    # traditional RF
    N = 300
    Xvals = np.arange(-15, 15, 10. / (2.0 * N))

    RF_DOG = []
    FFT_RF = []
    RField = []

    RF_DOG = findSpacing(Xvals, cone_spacing)
    FFT_RF = Fourier(RF_DOG, N)
    length = np.floor(FFT_RF.shape[0] / 2.) + 1
    RField = normRField(freqs, Xvals[length:] * 60,
                             FFT_RF[length:])

    receptive_field = {
                            'length': length,
                            'dog': RF_DOG,
                            'coneResponse': {'fft': FFT_RF, },
                            'fft': RField,
                            'xvals': Xvals
                            }
    return receptive_field


def findSpacing(Xvals, cone_spacing=2.0):
    SD = 0.1
    dog = DoG(Xvals, SD, 5.0)
    dist = findDist(Xvals, dog)
    i = 0

    print abs(cone_spacing - dist)
    if dist < cone_spacing:
        while dist < cone_spacing:
            SD += 0.05
            dog = DoG(Xvals, SD, 5.0)
            dist = findDist(Xvals, dog)
            i += 1
            if i == 99:
                raise IOError('Sorry, cannot find DoG')
    if dist > cone_spacing:
        while dist < cone_spacing:
            SD -= 0.01
            dog = DoG(Xvals, SD, 5.0)
            dist = findDist(Xvals, dog)
            i += 1
            if i == 99:
                raise IOError('Sorry, cannot find DoG')
    return dog


def findDist(Xvals, dog):

    mins = (np.diff(np.sign(np.diff(dog))) > 0).nonzero()[0] + 1
    print Xvals[mins[1]], Xvals[mins[0]]
    dist = Xvals[mins[1]] - Xvals[mins[0]]
    return dist


def sineWave(thisFrequency, xvals):
    """Generate a sine wave

    .. math::
       p = \\frac{1+ \\sin{(x*\lambda + \\frac{\\pi}{2})}}{2}


    with :math:`x` representing an array of locations in space (in the cone \
    receptive field), :math:`\\lambda` the spatial frequency, \
    :math:`\\frac{\\pi}{2}` ensures that the sine waves are in phase with the \
    receptive field and the remaining terms normalize the sine wave and bound \
    it between [0, 1].

    """
    # convert from cpd     (radians    / arcmin)
    Converted_freq = thisFrequency * (2 * np.pi) / 60.0
    sWave = (1.0 + np.sin(xvals * Converted_freq + (np.pi / 2.0))) / 2.0

    return sWave


def normRField(freqs, xp, yp):
    """return a normalized receptive field generated from a linear
    interpolation handle.
    """
    foo = np.interp(freqs[1:], xp, yp)
    rfield = foo / np.sum(foo)
    return rfield


def Fourier(recField, N):
    """return a normalized Fourier transformed receptive field.
    """
    FFT_RF = np.fft.fftshift(np.abs(np.fft.fft(recField))) / np.sqrt(2 * N)
    normFFT = FFT_RF / np.sum(FFT_RF)
    return normFFT


def DoG(xvals, excite_SD, inhibit_SD):
    """Generate a differenc of gaussian receptive field model.

    .. math::
       r(x) = \\frac{ \\exp{(-x^2)}}{2*excitatory_{SD}^2} -
       \\frac{ \\exp{(-x^2)}}{2*inhibitory_{SD}^2}


    where :math:`x` represents locations on the retina relative to a center \
    cone, and (:math:`excitatory_{SD}`) and (:math:`inhibitory_{SD}`) \
    represent the standard deviation of the excitatory center and inhibitory \
    surround, currently taken to be 0.5 and 5.0, respectively.

    """
    y_excite = gauss(xvals, excite_SD)
    y_inhibit = gauss(xvals, inhibit_SD)
    normFact = np.sum(y_excite) / np.sum(y_inhibit)
    y_inhibit *= normFact
    DoG_foo = y_excite - y_inhibit
    return DoG_foo / max(DoG_foo)


def gauss(x, SD):
    """A simple gaussian function
    """
    return 1.0 * np.exp(-(x) ** 2 / (2 * SD ** 2))


def normJayRField(freqs, jayRF):
    """return a normalized receptive field using a spline handle generated 
    using Jay's manual decomposition method.
    """
    foo = interpolate.splev(freqs[1:], jayRF, der = 0)
    rfield = foo / np.sum(foo)
    return rfield
