from __future__ import division
import matplotlib.pylab as plt
import sys
from eye.eyeModel import SchematicEye

from renderer import PlottingFun as pf

class Plotter(SchematicEye):
    """A plotting repo
    """
    def __init__(self, Analysis, recepitive_field, imageData, eyeOptics,
                 plots, save_plots, legend):
        # import our dictionaries of data:
        self.Analysis = Analysis
        self.rec_field = recepitive_field
        self.imageData = imageData
        self.eyeOptics = eyeOptics
        
        self.freqs = self.eyeOptics['freqs']

        #options:
        self.RFselect = self.rec_field['selection']        
        self.location = ['periph']
        if self.rec_field['fovea'] == True:
            self.location.append('fovea')
        # find where to save the plots based on sys architecture:
        self.findFigDirectory()
        
        # plot the appropriate plots, with options:
        if 'amp' in plots:
            self.plotAmpSpec(save_plots)
        if 'mtf' in plots:
            self.plotMTFfamily(save_plots)
        if 'accomm' in plots:
            self.plotMTFfamily(plot_option=2, save_plots=save_plots)
        if 'periph' in plots:
            self.plotPeripheral(save_plots)
        if 'plotDoG' in plots:
            self.plotDoG(save_plots)
        if 'plotDeconstructed' in plots:
            self.plotDeconstructedRF(save_plots)        
        if 'activity' in plots:
            self.plotActivity(save_plots)
        if 'info' in plots:
            self.plotInformation(save_plots, legend)
        
        

    def Find_PowerLaw_Placement(self):
        """Find where the :math:`\\frac{1}{f}` power law should be placed.
        
        .. note:: not yet working. Should find len(array)/2. Not much else.
        
        .. todo::
           1 / f heuristic
           
        """
        pass

    
    def findFigDirectory(self):
        """
        """

        if sys.platform == 'i686':
            self.figPath = 'C:/Users/Brian/Documents/eschaton/Figures/'
        if sys.platform == 'darwin':
            self.figPath = '~/Documents/GDrive/eschaton/Figures/'
    
    def plotInformation(self, save_plots=False, legend=False):
        """
        :param cones: list of cones in analysis
        :type cones: list
        :param save_plots: decide whether to save plots (True) or not (False)
        :type save_plots: bool
        :param legend: decide whether to add a legend (True) or not (False)
        :type legend: bool
        
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
            
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
    
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
        
        for key in self.Analysis:
            plt.plot(self.Analysis[key]['cones'], self.Analysis[key]['info'],
                     self.Analysis[key]['line']+'o', label=key,
                     linewidth=2.5, markersize=9)
    
    
        plt.xlim([min(self.Analysis[key]['cones'])-0.1, 
                  max(self.Analysis[key]['cones'])+0.1])
        if legend:
            plt.legend(loc='upper left')
        plt.xlabel('cones in mosaic')
        plt.ylabel('entropy (bits)')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'InfoPlot.png')
            plt.close()
        else:
            plt.show()
    
    def plotDoG(self, rec_field, save_plots=False):
    
        """
        Currently this function outputs the following:
    
        :param save_plots: decide whether to save plots (True) or not (False)
        :type save_plots: bool
    
        :returns: * Plot1: difference of gaussians receptive field \n
                  * Plot2: FFT spectrum of plot 1.
                
    
            
        .. figure:: ../../Figures/ConeRF.png 
           :height: 300px
           :width: 400px
           :align: center   
           
           **Fig 1:** A simple Diff of Gaussian receptive field
    
           
        .. figure:: ../../Figures/ConeRF_FFT.png  
           :height: 300px
           :width: 400px
           :align: center
           
           **Fig 2:** And converted into a probability density and Fourier 
           transformed
        
        """
        
        for loc in self.location:
            if self.RFselect == 'fft':
                fig = plt.figure(figsize=(8,6))
                ax = fig.add_subplot(111)
            
                pf.AxisFormat()
                pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
                
                ax.plot(self.rec_field['xvals'], self.rec_field['dog'][loc], 
                        'k', linewidth=2.5)
            
                plt.xlabel('distance (arcmin)')
                plt.ylabel('amplitude')
                
                plt.tight_layout()
                
                if save_plots:
                    fig.show()
                    fig.savefig(self.figPath + 'ConeRF' + loc + '.png')
                    plt.close()
                else:
                    plt.show()
                
                fig2 = plt.figure(figsize=(8,6))
                ax = fig2.add_subplot(111)
                pf.AxisFormat()
                pf.TufteAxis(ax, ['left', 'bottom'])
                
                length = self.rec_field['length']
                
                ax.loglog(self.rec_field['xvals'][length:] * 60, 
                    self.rec_field['coneResponse'][self.RFselect][loc][length:], 
                            'k', linewidth=2.5)
                
            
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
            
                plt.ylim([10**-4, 10**-1])
                plt.xlim([self.freqs[1], 100])
                plt.xlabel('spatial frequency (cycles / deg)')
                plt.ylabel('density')
                
                plt.tight_layout()
                
                if save_plots:
                    fig2.show()
                    fig2.savefig(self.figPath + 'ConeRF_FFT' + loc + '.png')
                    plt.close()
                else:
                    plt.show()
    
            
    def plotDeconstructedRF(self, save_plots):
        """
        :returns: two plots.
        
        .. figure:: ../../Figures/RecField_JAY.png
           :height: 300px
           :width: 400px
           :align: center
           
           **Fig 1:** Receptive field (black) with three spatial frequencies \
           (blue = 1 cpd, green = 5 cpd and red = 10 cpd)
    
    
        .. figure:: ../../Figures/FFT_JAY.png
           :height: 300px
           :width: 400px
           :align: center
           
           **Fig 2:** Response of receptive field to sine waves. 
           
        """
        for loc in self.location:    
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'], Nticks=[7,5])
            
            ax.plot(self.rec_field['xvals'], self.rec_field['sineWave'][:,1], 
                    'b-', linewidth =2, label='1 cpd')
            ax.plot(self.rec_field['xvals'], self.rec_field['sineWave'][:,5], 
                    'g-', linewidth = 2, label='5 cpd')
            ax.plot(self.rec_field['xvals'], self.rec_field['sineWave'][:,10], 
                    'r-', linewidth = 2, label='10 cpd')
                    
            ax.plot(self.rec_field['xvals'], self.rec_field['dog'][loc],
                    'k-', linewidth = 3)
            
            plt.xlim([min(self.rec_field['xvals']),
                      max(self.rec_field['xvals'])])
            plt.ylim([-0.15, 1.05])
        
            plt.xlabel('distance (arcmin)')
            plt.ylabel('amplitude')      
            
            plt.tight_layout()
        
            if save_plots:
                fig.show()
                fig.savefig(self.figPath + 'RecField_JAY' + loc + '.png')
                plt.close()
            else:
                plt.show()
        
            
            if save_plots:
                fig2 = plt.figure(figsize=(8,6))
                ax = fig2.add_subplot(111)
            else:
                fig = plt.figure(figsize=(8,6))
                ax = fig.add_subplot(111)
                
            pf.AxisFormat()
            pf.TufteAxis(ax, ['left', 'bottom'])
            
            length = self.rec_field['length']
            ax.loglog(self.rec_field['xvals'][length:] * 60,
                      self.rec_field['coneResponse']['jay'][loc][:], 
                        'k-', linewidth = 2.5)
            
            plt.ylim([10**-4, 10**-1])
            plt.xlim([self.freqs[1], 100])
            
            plt.xlabel('spatial frequency (cycles / deg)')
            plt.ylabel('density')      
            
            plt.tight_layout()
            
        
            if save_plots:
                fig2.show()
                fig2.savefig(self.figPath + 'FFT_JAY' + loc + '.png')
                plt.close()
            else:
                plt.show()
            
            
    def plotAmpSpec(self, save_plots = False):
        """Plot the amplitude spec of image.
    
        :param save_plots: decide whether to save plots (True) or not (False)
        :type save_plots: bool        
        
        **This function produces:**
        
        .. figure:: ../../Figures/ampSpec.png
           :height: 300px
           :width: 400px
           :align: center  
    
           **Fig 1:** An amplitude spectrum of a series of images.
           
        """
        
        plaw = self.imageData['powerlaw'](self.imageData['imagexval'])
            
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'])
    
        ax.loglog(self.imageData['imagexval'], self.imageData['ampMean'],
                    'r.', markersize = 10, markeredgewidth = 0)
     
        ax.loglog(self.imageData['imagexval'], plaw, 
                  'k', linewidth = 2.5)
        
        ax.text(1, 10**-2.4, r'$\frac{1}{\mathit{f}}$', size = 35)
    
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('density')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'ampSpec.png')
            plt.close()
        else:
            plt.show()
    
    
    def plotMTFfamily(self, save_plots = False, plot_option = 1, 
                       legend = False, legend_fig=False):    
        """Plot a family of MTF curves derived from a schematic eye.
        
        :param save_plots: decide whether to save plots (True) or not (False).
        :type save_plots: bool
        :param plot_option: decide what to plot,
                        1 = on axis
                        2 = off axis
        :param legend: turn legend on (True) or off (False). Default = True. \n
        :type legend: bool
        
        :returns: Modulation transfer function plots for a family of curves.
        
        .. note:: 
           * Need to add options for what to plot. 
           * Legend option is not currently working properly.
        
        **This produces:**
        
        .. figure:: ../../Figures/MTFfamily.png
           :height: 300px
           :width: 400px
           :align: center        
           
           **Fig 1:** A family of MTF curves from OSLO data.
           
        """
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
    
        ''' on axis plots '''
        ax.plot(self.freqs, self.eyeOptics['onAxis']['inf'][:,4], 'k', 
                linewidth = 2.5, label='diffraction ')
        ax.plot(self.freqs, self.eyeOptics['onAxis']['inf'][:,2], 'r', 
                linewidth=2.5, label='infinity') 
        
        if plot_option == 1:
            ax.plot(self.freqs, self.eyeOptics['onAxis']['20ft'][:,2], 'g', 
                    linewidth=2.5, label='20ft')
            ax.plot(self.freqs, self.eyeOptics['onAxis']['1m'][:,2], 'b',
                    linewidth=2.5, label='1m')
            ax.plot(self.freqs[:20], self.eyeOptics['onAxis']['16in'][:20,2], 
                    'm', linewidth=2.5, label='16in')
    
        ''' plot2 : off axis plots '''
    
        if plot_option == 2:
            ax.plot(self.freqs, self.eyeOptics['offAxis']['inf'][:,2], 
                    'r--', linewidth=2.5)
            
            ax.plot(self.freqs[:60], 
                    self.eyeOptics['object']['16in16in'][:60,2],
                    'g--', linewidth=2.5, label = 'near focus, near obj')
            ax.plot(self.freqs[:20], self.eyeOptics['object']['16under'][:20,2], 
                    'c--', linewidth=2.5, label = 'underacc, far obj')
            ax.plot(self.freqs[:20],
                    self.eyeOptics['object']['16in20ft'][:20,2], 
                    'b--',linewidth=2.5, label = 'near focus, far obj')
    
        if plot_option ==3:
            ax.plot(self.freqs, self.eyeOptics['offAxis']['inf'][:,4],
                    'k--', linewidth=2.5)
            
            ax.plot(self.freqs, self.eyeOptics['offAxis']['20ft'][:,2], 'g--', 
                    linewidth=2.5)
            ax.plot(self.freqs, self.eyeOptics['offAxis']['1m'][:,2], 'b--', 
                    linewidth=2.5)
            ax.plot(self.freqs[:20], self.eyeOptics['offAxis']['16in'][:20,2],
                    'm--', linewidth=2.5)
                    
        if legend: 
            ax.legend(loc='upper right')#title='object, retina')
        
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('modulation')
        
        plt.tight_layout()
        
        if save_plots:
            
            fig.show()
            if plot_option == 1:
                fig.savefig(self.figPath + 'MTFfamilyOnAxis.png')
            if plot_option == 2:
                fig.savefig(self.figPath + 'MTFfamily.png')
            plt.close()
        else:
            plt.show()
            
        """ Legend figure """
        if legend_fig:
            if save_plots:
                fig2 = plt.figure(figsize=(8,6))
                ax = fig2.add_subplot(111)
            else:
                fig = plt.figure(figsize=(8,6))
                ax = fig.add_subplot(111)
                
            font = {'weight' : 'norm', 'size'  : 32}
            legend = {'frameon' : False}
            plt.rc('font', **font)
            
            
            ax.plot(0,0, 'k', linewidth=2.5, label='diffraction')
            ax.plot(0,0, 'r', linewidth=2.5, label='infinity') 
            ax.plot(0,0, 'g--', linewidth=2.5, 
                    label = 'near focus, near object')
            ax.plot(0,0, 'c--', linewidth=2.5, label = 'underacc, far object')
            ax.plot(0,0, 'b--', linewidth=2.5,
                    label = 'near focus, far object')
            ax.legend(loc = 'center')
            plt.axis('off')
            
            if save_plots:
                fig2.show()
                fig2.savefig(self.figPath + 'Legend.png')
                plt.close()
            else:
                plt.show()
    
    
    
    
    def plotPeripheral(self, save_plots=False, legend=False):
        """
        Plot peripheral MTF with a comparison to Navarro et al 1993 or 
        Williams et al. 1996
        
        :param save_plots: decide whether to save plots (True) or not (False).
        :type save_plots: bool
        :param legend: turn legend on (True) or off (False). Default = True. \n
        :type legend: bool
        
        Currently supports 0, 10, 20, 40 degrees eccentricity.
        
        **This produces:**
        
        .. figure:: ../../Figures/MTFperiphery.png
           :height: 300px
           :width: 400px
           :align: center        
           
           **Fig 1:** A family of MTF curves from experimental data (dotted) 
           and schematic eye.        
        """
        from eye.Optics import MTF
        
        Fovea = MTF(self.freqs[:], 0)
        TenDeg = MTF(self.freqs[:], 10)
        TwentyDeg = MTF(self.freqs[:],20)
        FourtyDeg = MTF(self.freqs[:],40)
        
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'], [5,5])
        
        #Navarro et al 1993 analytical func:
        ax.plot(self.freqs[:], Fovea, 'm--', label='fovea', linewidth=2.5)
        ax.plot(self.freqs[:], TenDeg, 'r--', label='10 deg', linewidth=2.5)
        ax.plot(self.freqs[:], TwentyDeg, 'g--', label='20 deg', linewidth=2.5)
        ax.plot(self.freqs[:], FourtyDeg, 'b--',label='40 deg', linewidth=2.5)
        
        #OSLO ray trace data:
        ax.plot(self.freqs[:], self.eyeOptics['onAxis']['inf'][:,2], 
                'm-', label='fovea', linewidth=2.5)
        ax.plot(self.freqs[:], self.eyeOptics['offAxis']['inf'][:,2],
                'r-', label='10 deg', linewidth=2.5)        
        ax.plot(self.freqs[:], self.eyeOptics['farPeriph']['20deg'][:,2], 
                'g-', label='20 deg', linewidth=2.5)
        ax.plot(self.freqs[:30], self.eyeOptics['farPeriph']['40deg'][:30,2], 
                'b-', label='40 deg', linewidth=2.5)
                
        if legend: 
            ax.legend(loc='lower left')#,title='object dist, retinal location')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()
        
        #plt.ylim([10**-2.5, 10**0])
        #plt.xlim([self.freqs[1], 100])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('modulation')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'MTFperiphery.png')
            plt.close()
        else:
            plt.show()
            
    def plotActivity(self, Receptive_Field = 'FFT', save_plots = False, 
                   legend = False):
        """Plot powerlaw amplitude spectrum after accounting for MTF and MTF \
        + receptive field
        
        :param Receptive_Field: Indicate which receptive field to use. \
                                Option are 'FFT' or 'Jay'
        :type Receptive_Field: str        
        :param legend: turn legend on (True) or off (False). Default = True.
        :type legend: bool
        :param save_plots: decide whether to save a plot (True) or not (False)
        :type save_plots: bool
        :returns: Plot 1: power law times Modulation Transfer Function for 
                  three conditions \n
                  Plot 2: same plot with difference of gaussians receptive 
                  field included as well.  This is the final plot in this \
                  series.
        :rtype: plt.plot

        .. note::
           Eventually would like to introduce a heuristic to determine the 
           location of the power law text that is plotted. Currently hard \
           coded.
           
        
        
        **This function plots estimates of the photoreceptor activity**
        
        .. figure:: ../../Figures/MTFfamilyMod.png
           :height: 300px
           :width: 400px
           :align: center    
           
           **Fig 1:** Amplitude * MTF
        
        .. figure:: ../../Figures/MTFfamilyModAmp.png
           :height: 300px
           :width: 400px
           :align: center  
           
           **Fig2:** Amplitude * MTF * Receptive field
           
        """
            
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)

        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'])

        for key in self.Analysis:
            ax.loglog(self.Analysis[key]['freqs'],
                      self.Analysis[key]['preCone'],
                      self.Analysis[key]['line'],
                      linewidth=2.5, label=key)
        
        if legend: 
            ax.legend(loc='lower left')#,title='object dist, retinal location')
        
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        
        mi, ma = plt.ylim()
        
        ax.text(15, 10**-4.0, r'$\frac{1}{\mathit{f}}$', size = 35)
        
        plt.ylim([10**-8, 10**-3])
        plt.xlim([self.freqs[1], 100])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('density')
        
        plt.tight_layout()
        
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'MTFfamilyMod.png')
            plt.close()
        else:
            plt.show()
        
        """ 
        plot 2 
        """
       
        if save_plots:
            fig2 = plt.figure(figsize=(8,6))
            ax = fig2.add_subplot(111)
        else:
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            
        pf.AxisFormat()
        pf.TufteAxis(ax, ['left', 'bottom'])

        for key in self.Analysis:
            ax.loglog(self.Analysis[key]['freqs'],
                      self.Analysis[key]['retina'],
                      self.Analysis[key]['line'],
                      linewidth=2.5, label=key)                
        
        if legend: 
            ax.legend(loc='upper right')
        
        mi, ma = plt.ylim()

        ax.text(12, 10**-5.0, r'$\frac{1}{\mathit{f}}$', size = 35)
        
        plt.ylim([10**-8, 10**-3.75])
        plt.xlim([self.freqs[1], 100])
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('density')
        
        plt.tight_layout()
       
        if save_plots:
            if Receptive_Field.lower() == 'jay':
                save_name = self.figPath + 'MTFfamilyModAmpJAY.png'
            else:
                save_name = self.figPath + 'MTFfamilyModAmp.png'
            
            fig2.show()
            fig2.savefig(save_name)     
            plt.close()
            
        else:
            plt.show()
