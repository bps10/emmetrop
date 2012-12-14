import numpy as np
from scipy import interpolate


class ConeReceptiveFields(object):
    """This class generates models of cone receptive fields.
    
    Currently there is not a real need for a class. However, it is expected
    that as this project grows, we will introduce additional models and this
    class will surve a more dynamic purpose.
    """
    def __init__(self, freqs, excite_SD = None, inhibit_SD = None,
                 location = None):
        """
        """
        if not excite_SD:
            excite_SD = [0.5, 0.1]
        if not inhibit_SD:
            inhibit_SD = [5.0, 5.0]
        if not location:
            location = ['periph', 'fovea']
            
        self.genReceptiveFields(freqs, excite_SD, inhibit_SD, location)
        
    def genReceptiveFields(self, freqs, excite_SD, inhibit_SD, location ):
        """Create a difference of gaussians receptive field and plot it.
    
        :param excite_SD: standard deviation parameter of excitatory gaussian.
        :param inhibit_SD: standard deviation parameter of inhibitory \
        surround gaussian. Default = 5.0.   
        
        :returns: RF_DOG, RF_SPLINE, FFT_RF
        
        
    
        """
        
        # traditional RF
        N = 400
        Xvals = np.arange(-15, 15, 10./ (2.0 * N) )
        
        gauss1 = lambda x, excite_SD : 1.0*np.exp(-(x)**2 / (2 * excite_SD**2)) 
        gauss2 = lambda x, inhibit_SD : 1.0*np.exp(-(x)**2 / (2 * inhibit_SD**2))
        
        RF_DOG = {}
        FFT_RF = {}
        RF_SPLINE = {}
        RField = {}
        
        for i, loc in enumerate(location):
            
            y_excite = gauss1(Xvals, excite_SD[i])     
            y_inhibit = gauss2(Xvals, inhibit_SD[i]) 
    
            normFact = np.sum(y_excite) / np.sum(y_inhibit)
            y_inhibit *= normFact
        
    
            DoG_foo = y_excite - y_inhibit
            RF_DOG[loc] = DoG_foo / max(DoG_foo)
        
            FFT_RF[loc] = (np.fft.fftshift(np.abs(np.fft.fft(RF_DOG[loc]))) 
                        / np.sqrt(2 * N))
    
            length = np.floor(FFT_RF[loc].shape[0] / 2.) + 1    
        
            ## set up for interpolation
            RF_SPLINE[loc] = interpolate.splrep(Xvals[length:] * 60,
                                            FFT_RF[loc][length:], s=0)
            
            FFT_RF[loc] = FFT_RF[loc] / np.sum(FFT_RF[loc])
            
            #peripheral RF
            RField_foo = interpolate.splev(freqs[1:], RF_SPLINE[loc], der = 0)
            RField[loc] = RField_foo / np.sum(RField_foo)
            
    
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
                    
                # convert from cpd     (radians    / arcmin)
                Converted_freq = thisFrequency * (2 * np.pi) / 60.0 
                    
                sine_wave[:,i] = ( 1.0 + np.sin(Xvals * Converted_freq +
                                     (np.pi / 2.0) ))/ 2.0
                cone_response[i] = np.sum(sine_wave[:,i] * RF_DOG[loc])
    
            #symmetric_CR= np.zeros((len(cone_response)*2.)+2)
            #symmetric_CR[:len(cone_response)] = cone_response
            #symmetric_CR[len(cone_response)] = 10 #dc 
            #symmetric_CR[len(cone_response)+2:] = cone_response[::-1]
            CR = cone_response / np.sum(cone_response)
            Jay_CR[loc] = CR
            Jay_RF[loc] = interpolate.splrep(Xvals[length:] * 60,
                                                    cone_response, s=0)
            #print Jay_CR[loc]
            #foveal RF
            Jay_RField_foo = interpolate.splev(freqs[1:], Jay_RF[loc], der = 0)
            Jay_RField[loc] = Jay_RField_foo / np.sum(Jay_RField_foo)
            
            
    
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