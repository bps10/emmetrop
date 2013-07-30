from __future__ import division
import matplotlib.pylab as plt
import numpy as np

from base import plot as pf


def plotInformation(Analysis, save_plots=False, figpath='Figures/', legend=False):
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
    
    for key in Analysis:
        plt.plot(Analysis[key]['cones'], Analysis[key]['info'],
                 Analysis[key]['line']['style'],
                  c = Analysis[key]['line']['color'], label=key,
                 linewidth=2.5, markersize=9)


    plt.xlim([min(Analysis[key]['cones'])-0.1, 
              max(Analysis[key]['cones'])+0.1])
    if legend:
        plt.legend(loc='upper left')
    plt.xlabel('cones in mosaic')
    plt.ylabel('entropy (bits)')
    
    plt.tight_layout()
    
    if save_plots:
        fig.show()
        fig.savefig(figPath + 'InfoPlot.png')
        plt.close()
    else:
        plt.show()