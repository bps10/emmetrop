from analysis.NeitzModel import SchematicEyeAnalysis


def main():
    """Runs main sequence of analysis and generates plots in default mode
    
    .. todo::
       * Get user input to control plotting options.
       
    """
    
    Eye = SchematicEyeAnalysis()
    
    
    #Eye.MTFfamilyPlots(plot_option=2, legend = 0)
    #Eye.PlotPowerSpec()
    #Eye.DiffOfGaussian()  
    
    Eye.DeconstructedFFT(plot_option = 0)
    Eye.AmpModPlot(save_plots = False, Receptive_Field = 'Jay')
    
    return Eye
    
if __name__ == "__main__":
    
    Eye = main()
