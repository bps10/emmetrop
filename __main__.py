#!/usr/bin/env python
from emmetrop.analysis.NeitzModel import NeitzModel, genMeta
from emmetrop.renderer import plotRepo as pr
from emmetrop.scene.Images import Images

from base import cones as cones

import argparse


def main(args):
    """Runs main sequence of analysis and generates plots. This is for 
    command line control.
    """
    
    if args.save:
        save_plots = True
    else:
        save_plots = False
        
    analysis_args = []
    
    if args.distance:
        analysis_args.append('dist')
    if args.focus:
        analysis_args.append('focus')
    if args.off_axis or args.verbose:
        analysis_args.append('off_axis')
    if args.wavelength or args.verbose:
        analysis_args.append('wavelength')

    # get meta data:
    _meta, cpd = genMeta()

    # get image data stuff:
    ImageData = Images().returnImageData()

    # receptive field stuff:        
    rec_field = cones.genReceptiveFields(cpd, 2)

    Analysis, diffract = NeitzModel(ImageData, rec_field, cpd, _meta, analysis_args)



    plot_args = []

    if args.mtf or args.verbose:
        plot_args.append('mtf')   
    
    if args.amp or args.verbose:
        plot_args.append('amp')
      
    if args.activity or args.verbose:
        plot_args.append('activity')

    if args.dog or args.verbose:
        plot_args.append('plotDoG')

    if args.info:
        plot_args.append('info')

    if args.series:
        plot_args.append('seriesPlot')

    # send to renderer module for plotting:
    pr.Plotter(analysis_args, diffract, Analysis, rec_field, 
                ImageData, plot_args, save_plots=save_plots, legend=False)
        
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
    parser.add_argument("-q", "--series", action="store_true",
                        help="display wireframe series plot.")   
                        
    parser.add_argument("-d", "--distance", action="store_true",
                        help="include fovea in analyses")
    parser.add_argument('-f', "--focus", action="store_true",
                        help="include object set in analyses")
    parser.add_argument('-o', "--off_axis", action="store_true",
                        help="include far periphery in analyses")
    parser.add_argument('-w', "--wavelength", action="store_true",
                        help="include wavelength in analyses")

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="run all analyses and plots")
    parser.add_argument("-s", "--save", action="store_true",
                        help="save all plots")
                                                
    parser.add_argument("-e", "--eyegrow", action="store_true",
                        help="plot predicted eye growth against age")
                        
    args = parser.parse_args()
    
    Eye = main(args)
