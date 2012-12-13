from analysis.NeitzModel import SchematicEyeAnalysis
import argparse


def main(args):
    """Runs main sequence of analysis and generates plots in default mode
    
    .. todo::
       * Get user input to control plotting options.
       
    """
    
    
        
    if args.fovea:
        inc_fovea = True
    else:
        inc_fovea = False

        
    if args.save:
        save_plots = True
    else:
        save_plots = False


    if args.rfield.lower() == 'jay':
        Receptive_Field = "Jay"
    elif args.rfield.lower() == 'fft':
        Receptive_Field = "FFT"
    else: 
        raise('Sorry, receptive field entered not supported: Jay or FFT')
        
        
    Eye = SchematicEyeAnalysis(RF_opt = Receptive_Field, plot_args='all', 
                               save_arg = save_plots, fovea = inc_fovea)
    

        

        
    if args.mtf or args.verbose:
        Eye.MTFfamilyPlots(save_plots, plot_option=2, legend = 0)
        
    if args.periph or args.verbose:
        Eye.PeripheralPlot(save_plots)        
        
    if args.dog or args.verbose:
        if Receptive_Field.lower() == 'fft':
            Eye.PlotDoG(save_plots)
        else:
            pass
        
    if args.jay or args.verbose:
        if Receptive_Field.lower() == 'jay':
            plot = True
        else:
            plot = False
        Eye.DeconstructedFFT(plot_option = plot, save_plots = save_plots)
        
    if args.amp or args.verbose:
        Eye.PlotPowerSpec(save_plots)
      
    if args.activity or args.verbose:
        Eye.AmpModPlot(Receptive_Field, save_plots)
        
    if args.eyegrow:
        from analysis import Eye_Grow as eg
        eg.main()
    
    
    return Eye


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-m", "--mtf", action="store_true", 
                        help="plot MTF family")
    parser.add_argument("-a", "--amp", action="store_true",
                        help="plot amplitude spectrum")
    parser.add_argument("-d", "--dog", action="store_true", 
                         help="display DoG receptive field")
    parser.add_argument("-j", "--jay", action="store_true",
                        help="display deconstructed FFT (Jay)")
    parser.add_argument("-y", "--activity", action="store_true",
                        help="display activity plots")
    parser.add_argument("-p", "--periph", action="store_true", 
                        help="display a square of a given number")
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="run all analyses and plots")
    parser.add_argument("-s", "--save", action="store_true",
                        help="save all plots")
    parser.add_argument("-r", "--rfield", type=str, default = 'FFT',
                        help="change receptive field used for analysis")
    parser.add_argument("-f", "--fovea", action="store_true",
                        help="include fovea in analyses")
    parser.add_argument("-e", "--eyegrow", action="store_true",
                        help="plot predicted eye growth against age")
                        
    args = parser.parse_args()
    
    Eye = main(args)
