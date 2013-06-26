from __future__ import division
import numpy as np

# emmetrop imports:
from emmetrop.cones.dogRFields import genReceptiveFields
from emmetrop.eye.eyeModel import SchematicEye
from emmetrop.scene.Images import Images
from emmetrop.scene.DataManip import rad2deg
from emmetrop.analysis import Information as info
from emmetrop.renderer import plotRepo as pr


class SchematicEyeAnalysis(object):
    """This class is designed to estimate the response of a linear 
    photoreceptor to an amplitude spectrum derived from natural images.
    The modulation transfer functions are computed in OSLO ray tracing
    software.
    
    **Legend for figures:**
    
    .. figure:: ../../Figures/Legend.png
       :height: 200px
       :width: 300px
       :align: center
       
       **Fig 1: Use this legend for all figures in this section**
       
    .. todo::
       * single graphic with all subfigures?
       * better plotting options.
    """
    
    def __init__(self, analysis_args, plot_args, save_arg):
        """
        
        """
        # get ray tracer data:
        self.EyeOptics = SchematicEye()

        # get meta data:
        self.genMeta()

        # get image data stuff:
        self.ImageData = Images().returnImageData()
        
        # receptive field stuff:        
        self.rec_field = genReceptiveFields(self.freqs, 2)
        
        # analysis stuff:
        self.NeitzModel(analysis_args)
        
        # send to renderer module for plotting:
        pr.Plotter(self.freqs, self.Analysis, self.rec_field, self.ImageData, 
                    plot_args, save_arg, legend = False)

    def genMeta(self):
        '''Just temporary.
        '''
        # get from ray tracer eventually
        self._meta = {}
        self._meta['retImg'] = 2 #np.max(self.xvals) # size of image in mm 
        self._meta['eye_length_'] = 24.2
        radians = np.tan(self._meta['retImg'] / self._meta['eye_length_'])
        self._meta['deg'] = rad2deg(radians)
        self._meta['mm/deg'] = self._meta['deg'] / self._meta['retImg']

        #self.xvals = np.arange(0.0005, 0.2, 0.0005)
        cycles = (np.arange(0, 399)) / 2
        self.freqs = cycles * self._meta['mm/deg'] 

    def NeitzModel(self, analysis_args):
        '''This function organizes the entire operation.
        A dictionary self.Analysis is created to reflect the 
        user options. All subsequent methods will use the keys
        of this dictionary.
        '''
        self.Analysis = {}

        if 'distance' in analysis_args:
            dist_range = 10 ** (np.arange(5, 23) / 3.0)
        else:
            dist_range = np.array([1e7])

        if 'focus' in analysis_args:
            focus_range = np.arange(0, 2, 0.1)
        else: 
            focus_range = np.array([0])

        if 'off_axis' in analysis_args:
            axis_range = np.arange(0, 21, 2)
        else:
            axis_range = np.array([0])

        if 'pupil_size' in analysis_args:
            pupil_range = np.arange(2, 9, 1)
        else:
            pupil_range = np.array([4])

        j = 0
        for dist in dist_range:
            for focus in focus_range:
                for axis in axis_range:
                    for pupil in pupil_range:

                        self.Analysis[j] =  {
                                            'dist': dist,
                                            'focus': focus,
                                            'off_axis': axis,
                                            'pupil_size': pupil, 
                                            'line': self.addLineStyle(),
                                            }
                        j += 1
        
        self.ComputeConeActivity()      
        self.TotalActivity()
        self.estimateInfo()

    def addLineStyle(self):
        """Add line style for use with plots
        """
        axis = {'diffract': 'k', 'inf': 'r', '16in20ft': 'b', 
                '16in16in': 'g', '16inunder':'c', '20deg':'m', '40deg':'m'}
        state = {'onAxis': '-', 'offAxis': '--', 'object': '--',
                 'farPeriph': '-.'}
            
        #line = axis[self.Analysis[key]['params'][1]]
        #line += state[self.Analysis[key]['params'][0]]
        line = 'k-'
        return line
        
    def ComputeConeActivity(self, brownian_motion=True):
        """Compute the estimated activity of a cone photoreceptor.
        
        :param Receptive_Field: a handle to the spline fitted receptive field.
        :type Receptive_Field: function handle.
        
        .. note:: 
           very preliminary. Completely hard coded.
        
        """
        Rec_Field = self.rec_field['fft']  
            
        self.ImageData['fitLaw'] = self.ImageData['powerlaw'](
                                                 self.freqs[1:])

        powerlaw = self.ImageData['fitLaw']
        if brownian_motion:
            from emmetrop.eye.movement import brownian_motion
            temp = np.arange(1, 80)
            spat = self.freqs[1:]
            movement_filter = brownian_motion(spat, temp)
            powerlaw *= movement_filter
            #powerlaw = powerlaw / sum(powerlaw)
            
        for key in self.Analysis:
            ind = [0, 100] 

            # generate MTFs for each condition:
            intensity = self.EyeOptics.traceEye(
                                    self.Analysis[key]['dist'], 
                                    self.Analysis[key]['off_axis'], 
                                    self.Analysis[key]['pupil_size'], 
                                    self.Analysis[key]['focus'])
            self.Analysis[key]['mtf'] = self.EyeOptics.genMTF(intensity)

            self.Analysis[key]['preCone'] = (powerlaw[ind[0]:ind[1]] * 
                self.Analysis[key]['mtf'][ind[0]:ind[1]])

            self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] *
                                        Rec_Field[ind[0]:ind[1]])
                  
    def TotalActivity(self, print_opt = False):
        """Compute the estimated activity in a photoreceptor.
        
        :param self: no params yet.
        
        :returns: percent of integrated.
        
        .. note:: 
           Completely hard coded right now.  Would like to change that.
        
        """
        
        for key in self.Analysis:
            self.Analysis[key]['total'] = np.sum(self.Analysis[key]['retina'])
            #self.Analysis[key]['percent'] = (self.Analysis[key]['total'] /
            #                        self.Analysis['diffract periph']['total'])

        if print_opt:
            print ' '
            print 'Total activity (proportion of diffraction)'
            print '-------------------------------------------'
            for key in self.Analysis:
                print key, ': ', self.Analysis[key]['percent']


        
    def estimateInfo(self, print_opt=False):
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

        
        total_images = self.ImageData['totalImages']        

        for key in self.Analysis:
            self.Analysis[key]['cones'] = [1,2,3,4,5,6]
            self.Analysis[key]['info'] = np.zeros(len(
                                                self.Analysis[key]['cones']))

        for amp in self.ImageData['rawAmp']:

            for key in self.Analysis:
                ind = 100
                fooInfo = info.SingleConeEntropyFunc((amp[:ind]**2 *
                                            self.Analysis[key]['retina']), 
                                            self.Analysis[key]['cones'])      
                self.Analysis[key]['info'] += fooInfo / total_images


        if print_opt == True:
            print ' '
            print 'Information'
            print '------------'
            for key in self.Analysis:
                print key, ': ', self.Analysis[key]['info']


if __name__ == "__main__":
    eye = SchematicEyeAnalysis()