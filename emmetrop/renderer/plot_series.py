from __future__ import division
import matplotlib.pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

from base import plot as pf
from emmetrop.scene import SignalProcessing as sig

def plotSeries(cpd, Analysis, analysis_args, figPath='Figures/', save_plots=False):
    '''
    '''

    for analysis in analysis_args:
        ylabel = (analysis.replace('dist', '$log_{10}$ (mm)')
                    .replace('off_axis', 'angle (deg)')
                    .replace('pupil_size', 'pupil size (mm)')
                    .replace('focus', 'lens (D)'))

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')     
        pf.AxisFormat()

        keys = 0
        for key in Analysis: keys += 1
        samples = len(Analysis[0]['mtf'])
        X = np.zeros((50, keys))
        Y = np.zeros((50, keys))
        Z = np.zeros((50, keys))

        for key in Analysis:
            X[:, key] = cpd[:50]
            if analysis == 'dist':
                Y[:, key] = np.ones(50) * np.log10(Analysis[key]['dist'])
            else:
                Y[:, key] = np.ones(50) * Analysis[key][analysis]
            Z[:, key] = Analysis[key]['mtf'][:50]

        ax.plot_wireframe(X, Y, Z)
        ax.set_xlabel('cycles/degree')
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'MTFfamily.png')
            plt.close()
        else:
            plt.show()

        ## Retina plot ##
        #################

        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111, projection='3d')     
        pf.AxisFormat()

        for key in Analysis:
            X[:, key] = cpd[:50]
            if analysis == 'dist':
                Y[:, key] = np.ones(50) * np.log10(Analysis[key]['dist'])
            else:
                Y[:, key] = np.ones(50) * Analysis[key][analysis]
            contrast = (Analysis[key]['retina'][:50] / 
                    np.max(Analysis[key]['retina'][:50]))
            Z[:, key] = sig.decibels(contrast)

        ax.plot_wireframe(X, Y, Z)
        ax.set_xlabel('cycles/degree')
        ax.set_ylabel(ylabel)
        plt.tight_layout()
        if save_plots:
            fig.show()
            fig.savefig(self.figPath + 'MTFfamily.png')
            plt.close()
        else:
            plt.show()
