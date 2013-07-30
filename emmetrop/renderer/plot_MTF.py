from __future__ import division
import matplotlib.pylab as plt
import numpy as np

from base import plot as pf

def plotMTFfamily(cpd, Analysis, diffract, figPath='Figures/', 
	save_plots=False, legend=False):    
    """Plot a family of MTF curves derived from a schematic eye.
    
    :param save_plots: decide whether to save plots (True) or not (False).
    :type save_plots: bool
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
    
    for key in Analysis:
        mtf = Analysis[key]['mtf']
        ax.plot(cpd, mtf,
                  Analysis[key]['line']['style'],
                  c = Analysis[key]['line']['color'],)
    ax.plot(cpd, diffract['mtf'], 'k-', linewidth=2)

    if legend: 
        ax.legend(loc='upper right')#title='object, retina')

    #plt.ylim([self.min_dB, 1])
    plt.xlim([0, 60])
        
    plt.xlabel('spatial frequency (cycles / deg)')
    plt.ylabel('modulation transfer')
    
    plt.tight_layout()
    
    if save_plots:
        fig.show()
        fig.savefig(figPath + 'MTFfamily.png')
        plt.close()
    else:
        plt.show()