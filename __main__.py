#!/usr/bin/env python
from emmetrop.analysis.NeitzModel import NeitzModel, genMeta
from emmetrop import renderer as pr
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

    if analysis_args != []:
        # get meta data:
        _meta, cpd = genMeta()

        # get image data stuff:
        imageData = Images().returnImageData()

        # receptive field stuff:        
        rec_field = cones.genReceptiveFields(cpd, 2)

        Analysis, diffract = NeitzModel(imageData, rec_field, cpd, _meta, analysis_args)


        if args.mtf or args.verbose:
            pr.plotMTFfamily(cpd, Analysis, diffract, figPath='Figures/', 
                save_plots=False, legend=False)   
        
        if args.amp or args.verbose:
            pr.plotAmpSpec(imageData, figPath='Figures/', _brownian=True, 
                save_plots=save_plots)
          
        if args.activity or args.verbose:
            pr.plotActivity(cpd, Analysis, diffract, figPath='Figures/', 
                save_plots=save_plots, legend=False)

        if args.info:
            pr.plotInformation(Analysis, figpath='Figures/', 
                save_plots=save_plots, legend=False)

        if args.series:
            pr.plotSeries(cpd, Analysis, analysis_args, save_plots=False)

    if args.dog or args.verbose:
        _meta, cpd = genMeta()
        rec_field = cones.genReceptiveFields(cpd, 2)
        
        pr.plotDoG(rec_field, min_dB=20, figPath='Figures', save_plots=False)

    if args.big:
        pr.big_analysis_plot()    

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
                                              
    parser.add_argument("-b", "--big", action="store_true",
                        help="plot big analysis")                                              
    parser.add_argument("-e", "--eyegrow", action="store_true",
                        help="plot predicted eye growth against age")
                        
    args = parser.parse_args()
    
    Eye = main(args)
