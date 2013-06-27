from emmetrop.analysis.NeitzModel import SchematicEyeAnalysis
import argparse


def main(args):
    """Runs main sequence of analysis and generates plots. This is for 
    command line control.
    """
    
    if args.save:
        save_plots = True
    else:
        save_plots = False
        

    plot_args = []

    if args.mtf or args.verbose:
        plot_args.append('mtf')   
    
    if args.amp or args.verbose:
        plot_args.append('amp')
      
    if args.activity or args.verbose:
        plot_args.append('activity')
    
    if args.info or args.verbose:
        plot_args.append('info')

    if args.comp or args.verbose:
        plot_args.append('comp')

    if args.dog or args.verbose:
        plot_args.append('plotDoG')
    
    analysis_args = []
    
    if args.distance:
        analysis_args.append('distance')
    if args.focus:
        analysis_args.append('focus')
    if args.off_axis or args.verbose:
        plot_args.append('seriesPlot')
        analysis_args.append('off_axis')
    
    if not args.eyegrow:    
        Eye = SchematicEyeAnalysis(analysis_args=analysis_args,
                                   plot_args=plot_args, 
                                   save_arg=save_plots)
        return Eye
        
    if args.eyegrow:
        from emmetrop.analysis import Eye_Grow as eg
        eg.main()

    
    


if __name__ == "__main__": 
    parser = argparse.ArgumentParser()
    
    parser.add_argument("-m", "--mtf", action="store_true", 
                        help="display MTF family plots")
    parser.add_argument("-a", "--amp", action="store_true",
                        help="display amplitude spectrum plot")
    parser.add_argument("-r", "--dog", action="store_true", 
                         help="display DoG receptive field")
    parser.add_argument("-y", "--activity", action="store_true",
                        help="display activity plots")
    parser.add_argument("-i", "--info", action="store_true",
                        help="display information plot")
    parser.add_argument("-c", "--comp", action="store_true",
                        help="display comparison plot with Williams et al.")
                        
    parser.add_argument("-d", "--distance", action="store_true",
                        help="include fovea in analyses")
    parser.add_argument('-f', "--focus", action="store_true",
                        help="include object set in analyses")
    parser.add_argument('-o', "--off_axis", action="store_true",
                        help="include far periphery in analyses")
                        
    parser.add_argument("-v", "--verbose", action="store_true",
                        help="run all analyses and plots")
    parser.add_argument("-s", "--save", action="store_true",
                        help="save all plots")
                                                
    parser.add_argument("-e", "--eyegrow", action="store_true",
                        help="plot predicted eye growth against age")
                        
    args = parser.parse_args()
    
    Eye = main(args)
