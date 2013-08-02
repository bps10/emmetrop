from __future__ import division
import matplotlib.pylab as plt
import numpy as np

from base import plot as pf
from emmetrop.scene import SignalProcessing as sig
            
def plotActivity(cpd, Analysis, diffract, figPath='Figures/', 
        save_plots=False, legend=False):
    """Plot powerlaw amplitude spectrum after accounting for MTF and MTF \
    + receptive field
           
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
    pf.TufteAxis(ax, ['left', 'bottom'], [5, 5])

    for key in Analysis:
        contrast = (Analysis[key]['preCone'] / 
                np.max(diffract['preCone']))

        ax.semilogx(cpd[0:100], sig.decibels(contrast),
                  Analysis[key]['line']['style'],
                  c = Analysis[key]['line']['color'], label=key)
    contrast = diffract['preCone'] / np.max(diffract['preCone'])
    ax.plot(cpd[0:100], sig.decibels(contrast),
        'k-', linewidth=2)

    if legend: 
        ax.legend(loc='lower left')
    
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()
   
    plt.ylim([-20, 0])
    plt.xlim([0, 60])
    
    plt.xlabel('spatial frequency (cycles / deg)')
    plt.ylabel('contrast sensitivity (dB)')
    
    plt.tight_layout()
    
    if save_plots:
        fig.show()
        fig.savefig(figPath + 'MTFfamilyMod.png')
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
    pf.TufteAxis(ax, ['left', 'bottom'], [5, 5])

    for key in Analysis:
        contrast = (Analysis[key]['retina'] / 
                np.max(diffract['retina']))
        ax.semilogx(cpd[0:100], sig.decibels(contrast),
                  Analysis[key]['line']['style'],
                  c = Analysis[key]['line']['color'], label=key)
    contrast = diffract['retina'] / np.max(diffract['retina'])                
    ax.plot(cpd[0:100], sig.decibels(contrast),
                'k-', linewidth=2)

    if legend: 
        ax.legend(loc='upper right')
    
    mi, ma = plt.ylim()

    plt.ylim([-20, 0])
    plt.xlim([0, 60])
    
    plt.xlabel('spatial frequency (cycles / deg)')
    plt.ylabel('contrast sensitivity (dB)')
    
    plt.tight_layout()
   
    if save_plots:
        save_name = figPath + 'MTFfamilyModAmp.png'
        
        fig2.show()
        fig2.savefig(save_name)     
        plt.close()
        
    else:
        plt.show()