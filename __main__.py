from eschaton.analysis.NeitzModel import SchematicEyeAnalysis
import argparse


def main(args):
    """Runs main sequence of analysis and generates plots. This is for 
    command line control.
    """
    
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
        

    plot_args = []

    if args.mtf or args.verbose:
        plot_args.append('mtf')
        
    if args.periphPlot or args.verbose:
        plot_args.append('periph')        
    
    if args.accomm or args.verbose:
        plot_args.append('accomm')
        
    if args.dog or args.verbose:
        if Receptive_Field.lower() == 'fft':
            plot_args.append('plotDoG')
        elif Receptive_Field.lower() == 'jay':
            plot_args.append('plotDeconstructed')
                
    if args.amp or args.verbose:
        plot_args.append('amp')
      
    if args.activity or args.verbose:
        plot_args.append('activity')
    
    if args.info or args.verbose:
        plot_args.append('info')
    
    analysis_args = []
    
    if args.fovea or args.verbose:
        analysis_args.append('fovea')
    if args.objectSet or args.verbose:
        analysis_args.append('objectSet')
    if args.farPeriph or args.verbose:
        analysis_args.append('farPeriph')
    
    if not args.eyegrow:    
        Eye = SchematicEyeAnalysis(RF_opt = Receptive_Field, 
                                   analysis_args = analysis_args,
                                   plot_args=plot_args, 
                                   save_arg = save_plots)
        return Eye
        
    if args.eyegrow:
        from eschaton.analysis import Eye_Grow as eg
        eg.main()

    
    


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-c", "--accomm", action="store_true",
                        help="display accommodation MTF plots")
    parser.add_argument("-m", "--mtf", action="store_true", 
                        help="display MTF family plots")
    parser.add_argument("-a", "--amp", action="store_true",
                        help="display amplitude spectrum plot")
    parser.add_argument("-d", "--dog", action="store_true", 
                         help="display DoG receptive field")
    parser.add_argument("-y", "--activity", action="store_true",
                        help="display activity plots")
    parser.add_argument("-t", "--periphPlot", action="store_true", 
                        help="display peripheral mtf plots")
    parser.add_argument("-i", "--info", action="store_true",
                        help="display information plot")
                        
    parser.add_argument("-f", "--fovea", action="store_true",
                        help="include fovea in analyses")
    parser.add_argument('-o', "--objectSet", action="store_true",
                        help="include object set in analyses")
    parser.add_argument('-p', "--farPeriph", action="store_true",
                        help="include far periphery in analyses")
                        
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="run all analyses and plots")
    parser.add_argument("-s", "--save", action="store_true",
                        help="save all plots")
                                                
    parser.add_argument("-e", "--eyegrow", action="store_true",
                        help="plot predicted eye growth against age")
    parser.add_argument("-r", "--rfield", type=str, default = 'FFT',
                        help="change receptive field used for analysis")
                        
    args = parser.parse_args()
    
    Eye = main(args)
