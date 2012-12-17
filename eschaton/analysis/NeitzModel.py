from __future__ import division
import numpy as np

# eschaton imports:
from eschaton.cones.dogRFields import ConeReceptiveFields
from eschaton.eye.eyeModel import SchematicEye
from eschaton.scene.Images import Images

from eschaton.analysis import Information as info
from eschaton.renderer import plotRepo as pr

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
    
    def __init__(self, RF_opt, analysis_args, plot_args, save_arg):
        """
        
        """
        
        # get ray tracer data:
        self.EyeOptics = SchematicEye().returnOSLOdata()
        
        # get image data stuff:
        self.ImageData = Images().returnImageData()
        
        # receptive field stuff:        
        self.rec_field = ConeReceptiveFields(
                            self.EyeOptics['freqs']).returnReceptiveField()
        
        # append user options:
        self.rec_field['selection'] = RF_opt.lower()
        
        # analysis stuff:
        self.NeitzModel(analysis_args, RF_opt)
        
        # send to renderer module for plotting:
        pr.Plotter(self.Analysis, self.rec_field, self.ImageData, 
                   self.EyeOptics, plot_args, save_arg, legend = False)

    def NeitzModel(self, analysis_args, RF_opt):
        """This function organizes the entire operation...
        """
        self.Analysis = {}
        
        end = len(self.EyeOptics['offAxis']['diffract'])

        # Periph
        self.Analysis['diffract periph'] = {'ind': [1,end,60],
                                'params': ['offAxis','diffract', 'periph']}
        self.Analysis['inf perip'] = {'ind': [1,end,60],
                                'params': ['offAxis','inf', 'periph']}

        if 'fovea' in analysis_args:        
            self.Analysis['diffract fovea'] = {'ind': [1,end,60],
                                'params': ['onAxis','diffract', 'fovea']}
            self.Analysis['inf fovea'] = {'ind': [1,end,60],
                                'params': ['onAxis','inf', 'fovea']}  
                                
        if 'objectSet' in analysis_args:                        
            self.Analysis['near focus, far object'] = {'ind': [1,8,7],
                                    'params': ['object','16in20ft', 'periph']}
            self.Analysis['near focus, near object'] = {'ind': [1,60,59],
                                    'params': ['object','16in16in', 'periph']}
            self.Analysis['underaccomm, far object'] = {'ind': [1,8,7],
                                    'params': ['object','16inunder', 'periph']}
    
        if 'farPeriph' in analysis_args:
            self.Analysis['20deg periph'] = {'ind': [1,end,60],
                                'params': ['farPeriph', '20deg', 'periph']}
            self.Analysis['40deg periph'] = {'ind': [1,end,60],
                                'params': ['farPeriph', '40deg', 'periph']}

        self.addLineStyle()
        # run it all:
        self.ComputeConeActivity(RF_opt)       
        self.TotalActivity()
        self.estimateInfo(RF_opt)        

    def addLineStyle(self):
        """Add line style for use with plots
        """
        axis = {'diffract': 'k', 'inf': 'r', '16in20ft': 'b', 
                '16in16in': 'g', '16inunder':'c', '20deg':'m', '40deg':'m'}
        state = {'onAxis': '-', 'offAxis': '--', 'object': '--',
                 'farPeriph': '-.'}
        
        for key in self.Analysis:
            
            line = axis[self.Analysis[key]['params'][1]]
            line += state[self.Analysis[key]['params'][0]]
            self.Analysis[key]['line'] = line
        
        
    def ComputeConeActivity(self, Receptive_Field):
        """Compute the estimated activity of a cone photoreceptor.
        
        :param Receptive_Field: a handle to the spline fitted receptive field.
        :type Receptive_Field: function handle.
        
        .. note:: 
           very preliminary. Completely hard coded.
        
        """
        if Receptive_Field.lower() == 'jay':
            Rec_Field = self.rec_field['RField']['jay']
            print 'Used Jay receptive field'
            
        elif Receptive_Field.lower() == 'fft':
            Rec_Field = self.rec_field['RField']['fft']
            print 'Used FFT receptive field'
            
        else :
            Rec_Field = self.rec_field['RField']['ftt']
            
            print 'receptive field not understood, see options. Using FFT.'            
            
        self.ImageData['fitLaw'] = self.ImageData['powerlaw'](
                                                 self.EyeOptics['freqs'][1:])
        powerlaw = self.ImageData['fitLaw']
        
        for key in self.Analysis:
            ind = self.Analysis[key]['ind']
            LocKey = self.Analysis[key]['params'][0]
            OptKey = self.Analysis[key]['params'][1]
            RfKey = self.Analysis[key]['params'][2]
            
            self.Analysis[key]['preCone'] = (powerlaw[ind[0]-1:ind[1]-1]*
                                self.EyeOptics[LocKey][OptKey][ind[0]:ind[1]])
            self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] *
                                        Rec_Field[RfKey][ind[0]-1:ind[1]-1])
            self.Analysis[key]['freqs'] =self.EyeOptics['freqs'][ind[0]:ind[1]]
            
   

        


                  
    def TotalActivity(self, print_opt = True):
        """Compute the estimated activity in a photoreceptor.
        
        :param self: no params yet.
        
        :returns: percent of integrated.
        
        .. note:: 
           Completely hard coded right now.  Would like to change that.
        
        """
        
        for key in self.Analysis:
            self.Analysis[key]['total'] = np.sum(self.Analysis[key]['retina'])
            self.Analysis[key]['percent'] = (self.Analysis[key]['total'] /
                                    self.Analysis['diffract periph']['total'])

        if print_opt:
            print ' '
            print 'Total activity (proportion of diffraction)'
            print '-------------------------------------------'
            for key in self.Analysis:
                print key, ': ', self.Analysis[key]['percent']


        
    def estimateInfo(self, Receptive_Field, print_opt=False):
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
                ind = self.Analysis[key]['ind'][2]
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