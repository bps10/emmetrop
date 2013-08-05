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
    # get meta data:
    _meta, cpd = genMeta()

    if args.save:
        save_plots = True
    else:
        save_plots = False
        
    analysis_args = []
    rec_field = {}
    
    if args.distance:
        analysis_args.append('dist')
    if args.focus:
        analysis_args.append('focus')
    if args.off_axis or args.verbose:
        analysis_args.append('off_axis')
    if args.wavelength or args.verbose:
        _max = 0
        for wv in range(430, 701, 10):
            rec_field[wv] = cones.genReceptiveFields(cpd, 2, [559, 530], wv)

            length = rec_field[wv]['length']
            rec_sum = max(rec_field[wv]['coneResponse'][length:])
            if rec_sum > _max:
                _max = rec_sum
        rec_field['max'] = _max
        analysis_args.append('wavelength')

    else:
        rec_field[540] = cones.genReceptiveFields(cpd, 0.5)
        length = rec_field[540]['length']
        rec_field['max'] = max(rec_field[540]['coneResponse'][length:])

    if analysis_args != []:

        # get image data stuff:
        imageData = Images().returnImageData()

        Analysis, diffract = NeitzModel(imageData, rec_field, cpd, 
            args.field_angle, _meta, analysis_args, glasses=args.glasses)


        if args.mtf or args.verbose:
            pr.plotMTFfamily(diffract['cpd'], Analysis, diffract, figPath='Figures/', 
                save_plots=False, legend=False)   
        
        if args.amp or args.verbose:
            pr.plotAmpSpec(imageData, figPath='Figures/', _brownian=True, 
                save_plots=save_plots)
          
        if args.activity or args.verbose:
            pr.plotActivity(diffract['cpd'], Analysis, diffract, figPath='Figures/', 
                save_plots=save_plots, legend=False)

        if args.info:
            pr.plotInformation(Analysis, figpath='Figures/', 
                save_plots=save_plots, legend=False)

        if args.series:
            pr.plotSeries(diffract['cpd'], diffract, Analysis, 
                analysis_args, save_plots=False)

    if args.dog or args.verbose:

        rec_field = {}
        if args.wavelength:
            _max = 0
            for wv in range(430, 701, 10):
                rec_field[wv] = cones.genReceptiveFields(cpd, 2, [559, 530], wv)

                length = rec_field[wv]['length']
                rec_sum = max(rec_field[wv]['coneResponse'][length:])
                if rec_sum > _max:
                    _max = rec_sum
            rec_field['max'] = _max

        else:
            rec_field[1] = cones.genReceptiveFields(cpd, 2)
            length = rec_field[1]['length']
            rec_field['max'] = max(rec_field[1]['coneResponse'][length:])
        
        pr.plotDoG(rec_field, min_dB=20, figPath='Figures', save_plots=False)

    if args.big:
        pr.big_analysis_plot(args.glasses)
        if args.glasses:
            pr.glasses_comp_plot()    

    if args.eyegrow:
        from emmetrop.analysis import Eye_Grow as eg
        eg.main()

    if args.filter:
        pr.plotMovement()
        pr.plotGlasses()

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

    parser.add_argument('-z', '--field_angle', default='10', type=float)

    parser.add_argument("-v", "--verbose", action="store_true",
                        help="run all analyses and plots")
    parser.add_argument("-s", "--save", action="store_true",
                        help="save all plots")
                                              
    parser.add_argument("-b", "--big", action="store_true",
                        help="plot big analysis")      
    parser.add_argument("--glasses_analysis", action="store_true",
                        help="use big plot to plot diff between states")                                        
    parser.add_argument("-e", "--eyegrow", action="store_true",
                        help="plot predicted eye growth against age")
    parser.add_argument("--glasses", action="store_true",
                        help="add diffuser glasses")
    parser.add_argument("--filter", action="store_true",
                        help="plot movement and glasses filter")

    args = parser.parse_args()
    
    Eye = main(args)
