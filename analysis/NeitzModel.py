from __future__ import division
import numpy as np

# eschaton imports:
from cones.dogRFields import ConeReceptiveFields
from eye.eyeModel import SchematicEye
from scene.Images import Images

from analysis import Information as info
from renderer import plotRepo as pr

class SchematicEyeAnalysis(object):
    """This class is designed to estimate the response of a linear 
    photoreceptor to an amplitude spectrum derived from natural images.
    The modulation transfer functions are computed in OSLO ray tracing software.
    
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
    
    def __init__(self, RF_opt, plot_args, save_arg, fovea = True):
        """
        
        """
        
        # get ray tracer data:
        self.EyeOptics = SchematicEye().returnOSLOdata()
        
        # get image data stuff:
        self.ImageData = Images().returnImageData()
        
        # receptive field stuff:        
        self.rec_field = ConeReceptiveFields(self.EyeOptics['freqs']).returnReceptiveField()
        
        # append user options:
        self.rec_field['selection'] = RF_opt.lower()
        self.rec_field['fovea'] = fovea
        
        # analysis stuff:
        self.NeitzModel(fovea, RF_opt)
        
        # send to renderer module for plotting:
        pr.Plotter(self.Analysis, self.rec_field, self.ImageData, 
                   self.EyeOptics, plot_args, save_arg, legend = False)

    def NeitzModel(self, fovea, RF_opt):
        """
        """
        self.Analysis = {}
        # move to a method. Create keys, line styles by appending other keys.
        self.Analysis['diffract periph'] = {'line':'k-'}
        self.Analysis['inf perip'] = {'line':'r-'}
        if fovea:        
            self.Analysis['diffract fovea'] = {'line':'k-.'}
            self.Analysis['inf fovea'] = {'line':'r-.'}            
        self.Analysis['near focus, far object'] = {'line':'b--'}
        self.Analysis['near focus, near object'] = {'line':'g--'}
        self.Analysis['underaccomm, far object'] = {'line':'c--'}

        # run it all:
        self.ComputeConeActivity(RF_opt, fovea)       
        self.TotalActivity()
        self.estimateInfo(RF_opt, fovea)        

        
    def ComputeConeActivity(self, Receptive_Field, fovea = False):
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
            
        self.ImageData['fitLaw'] = self.ImageData['powerlaw'](self.EyeOptics['freqs'][1:])
        powerlaw = self.ImageData['fitLaw']
        
        for key in self.Analysis:
            if key == 'diffract periph':
                self.Analysis[key]['preCone'] = powerlaw * self.EyeOptics['offAxis']['diffract'][1:]
                self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] * 
                                                Rec_Field['periph'])
                self.Analysis[key]['freqs'] = self.EyeOptics['freqs'][1:]

            if key == 'diffract fovea':
                self.Analysis[key]['preCone'] = powerlaw * self.EyeOptics['onAxis']['diffract'][1:]
                self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] *
                                                Rec_Field['fovea'])
                self.Analysis[key]['freqs'] = self.EyeOptics['freqs'][1:]
                
            if key == 'inf perip':
                self.Analysis[key]['preCone'] = (powerlaw * 
                                            self.EyeOptics['offAxis']['inf'][1:,2])
                self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] * 
                                                Rec_Field['periph'])
                self.Analysis[key]['freqs'] = self.EyeOptics['freqs'][1:]        

            if key == 'inf fovea':
                self.Analysis[key]['preCone'] = (powerlaw * 
                                            self.EyeOptics['onAxis']['inf'][1:,2])
                self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] * 
                                                Rec_Field['fovea'])
                self.Analysis[key]['freqs'] = self.EyeOptics['freqs'][1:]  
                 
            if key == 'near focus, near object':
                self.Analysis[key]['preCone'] = (powerlaw[:59] * 
                                 self.EyeOptics['object']['16in16in'][1:60,2])
                self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] * 
                                                Rec_Field['periph'][:59])
                self.Analysis[key]['freqs'] = self.EyeOptics['freqs'][1:60]
            
            if key == 'near focus, far object':
                self.Analysis[key]['preCone'] = (powerlaw[:7] * 
                                 self.EyeOptics['object']['16in20ft'][1:8,2])                
                self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] * 
                                                Rec_Field['periph'][:7])
                self.Analysis[key]['freqs'] = self.EyeOptics['freqs'][1:8]        
            
            if key == 'underaccomm, far object':
                self.Analysis[key]['preCone'] = (powerlaw[:7] * 
                                self.EyeOptics['object']['16under'][1:8,2]) 
                self.Analysis[key]['retina'] = (self.Analysis[key]['preCone'] * 
                                                Rec_Field['periph'][:7])
                self.Analysis[key]['freqs'] = self.EyeOptics['freqs'][1:8]        

        


                  
    def TotalActivity(self):
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


        print ' '
        print 'Total activity (proportion of diffraction)'
        print '-------------------------------------------'
        for key in self.Analysis:
            print key, ': ', self.Analysis[key]['percent']


        
    def estimateInfo(self, Receptive_Field, fovea = False, print_option=False):
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

        if Receptive_Field.lower() == 'jay':
            Rec_Field = self.rec_field['RField']['jay']
            
        elif Receptive_Field.lower() == 'fft':
            Rec_Field = self.rec_field['RField']['fft']
            
        else :
            Rec_Field = self.rec_field['RField']['ftt']
        
        
        total_images = self.ImageData['totalImages']        

        for key in self.Analysis:
            self.Analysis[key]['cones'] = [1,2,3,4,5,6]
            self.Analysis[key]['info'] = np.zeros(len(self.Analysis[key]['cones']))

        for amp in self.ImageData['rawAmp']:

            for key in self.Analysis:
                
                #Diffraction:
                if key == 'diffract periph':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 *
                                                self.Analysis[key]['retina'] * 
                                                Rec_Field['periph']), 
                                                self.Analysis[key]['cones'])
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'diffract fovea':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 *
                                                self.Analysis[key]['retina'] * 
                                                Rec_Field['fovea']), 
                                                self.Analysis[key]['cones'])
                    self.Analysis[key]['info'] += fooInfo / total_images
           
                #Infinity
                if key == 'inf perip':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 * 
                                                self.Analysis[key]['retina'] * 
                                                Rec_Field['periph']),
                                                self.Analysis[key]['cones'])
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'inf fovea':
                    fooInfo = info.SingleConeEntropyFunc((amp[:60]**2 * 
                                                self.Analysis[key]['retina'] * 
                                                    Rec_Field['fovea']), 
                                                    self.Analysis[key]['cones'])
                    self.Analysis[key]['info'] += fooInfo / total_images
    
                if key == 'near focus, far object':                           
                    fooInfo = info.SingleConeEntropyFunc((amp[1:8]**2 *
                                            self.Analysis[key]['retina'] * 
                                            Rec_Field['periph'][1:8]), 
                                            self.Analysis[key]['cones'])
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'near focus, near object':            
                    fooInfo = info.SingleConeEntropyFunc((amp[1:60]**2 * 
                                            self.Analysis[key]['retina'] * 
                                            Rec_Field['periph'][1:60]), 
                                            self.Analysis[key]['cones'])
                    self.Analysis[key]['info'] += fooInfo / total_images

                if key == 'underaccomm, far object':                          
                    fooInfo = info.SingleConeEntropyFunc((amp[1:8]**2 *
                                            self.Analysis[key]['retina'] * 
                                            Rec_Field['periph'][1:8]), 
                                            self.Analysis[key]['cones'])
                                            
                    self.Analysis[key]['info'] += fooInfo / total_images 

        if print_option == True:
            print ' '
            print 'Information'
            print '------------'
            for key in self.Analysis:
                print key, ': ', self.Analysis[key]['info']
            
            



if __name__ == "__main__":
    eye = SchematicEyeAnalysis()