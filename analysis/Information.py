from __future__ import division
import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt

from scene import DataManip as dm
from scene import SignalProcessing as sig
from analysis import PlottingFun as pf


def SignalGen(DIRECTORY='C:/Data/UPenn_Images/LMS_Images_Raw', 
        AUX_DIRECTORY='C:/Data/UPenn_Images/AUX_Image_Files', Analysis_Option = 1, BIN_RESOLUTION = 300):
    """
    SignalGen(DIRECTORY,AUX_DIRECTORY,Analysis_Option)
    
    Generate LMS signals for Info analysis from LMS isomerization rates from
    the specified directories.  Used AUX_DIRECTORY to avoid files that have
    greater than 5# saturated pixels (as described in Garrigan et al.).  
    

    :param DIRECTORY: Location of Images to be analyzed.
    :param AUX_DIRECTORY: Location of Auxilarly files with warning info.
    
    :param Analysis_Option: 1 = R* histogram, 
                     2 = R*, Power Spectrum, 
                     3 = R*, PS & Correlation. 
                     Default = 1.
    :param BIN_RESOLUTION: Number of bins to use in histogram. Default = 300.    

    :returns: dict containing analysis of L,M,S cone isomerization rates and bin edges used to compute histogram..


    .. note::  
       * Still not quite right. May need to add additional image selection
         heuristics or look into normalizing for variance.  
       * Only option 1 has been tested rigorously.
       * can undo the cone settings by dividing by the scalar values obtained for
       * cone parameters: Tkacik et al. pg. 10

    """



    # First find an index of images that contain pixel saturations:
    try:
        AUX_Files = dm.getAllFiles(AUX_DIRECTORY,subdirectories = 2)
    except ValueError:
        print "No Auxilary files found. Check directory path or download from: \
        ftp://tofu.psych.upenn.edu/"
    goodFiles = np.zeros((len(AUX_Files)))
    for i in range(0, len(AUX_Files)):
        AUX_File = sio.loadmat(AUX_Files[i])
        if AUX_File['Image']['warning'][0][0].shape[0] > 0:
            goodFiles[i] = 0
        else :
            goodFiles[i] = 1


    #set up memory and variables used in all cases
    try: 
        files = dm.getAllFiles(DIRECTORY, subdirectories = 2)
    except ValueError:
        print "No images found. Check directory path or download from:\
        ftp://tofu.psych.upenn.edu/"
        
    stop = len(files)
    j = 0

    BINS = np.linspace(0,4 * 10**4, BIN_RESOLUTION + 1 )
    L_count = np.zeros((np.sum(goodFiles), len(BINS) - 1))
    M_count = np.zeros((np.sum(goodFiles), len(BINS) - 1))
    S_count = np.zeros((np.sum(goodFiles), len(BINS) - 1))

    # find the desired analysis option:

    if Analysis_Option == 1:
        for i in range(0,stop):
            if goodFiles[i] == 1:
                
                image = sio.loadmat(files[i]) 

                L_raw = image['LMS_Image'][:,:,0]
                L_raw = L_raw.flatten() / 100.0 # change into units R*/10ms
             
                M_raw = image['LMS_Image'][:,:,1]
                M_raw = M_raw.flatten() / 100.0

                S_raw = image['LMS_Image'][:,:,2]
                S_raw = S_raw.flatten() / 100.0
                
                L_count[j,:],bins = np.histogram(L_raw,BINS) 
                M_count[j,:],bins = np.histogram(M_raw,BINS)
                S_count[j,:],bins = np.histogram(S_raw,BINS)
                
                j += 1

        out = [{}]*4
        out[0] =   {
                     'name':'L',
                     'count':L_count
                     }
        out[1] =   {
                     'name':'M',
                     'count':M_count
                     }
        out[2] =     {
                     'name':'S',
                     'count':S_count
                     }
        out[3] =     {
                    'bins': BINS
                    }
        
    if Analysis_Option == 2:
    
        SIZE = sio.loadmat(files[0])
        SIZE = SIZE['LMS_Image'][:,:,0]
        SIZE = min(SIZE.shape)
        L_power = np.zeros((np.sum(goodFiles), np.floor(SIZE / 2 + 1)))
        M_power = np.zeros((np.sum(goodFiles), np.floor(SIZE / 2 + 1)))
        S_power = np.zeros((np.sum(goodFiles), np.floor(SIZE / 2 + 1)))

        for i in range(0,stop):
            if goodFiles[i] == 1:
                
                image = sio.loadmat(files[i]) 

                L_raw = image['LMS_Image'][:,:,0]
                L_power[j,:] = sig.welch2d(L_raw[0:SIZE,0:SIZE])
                L_raw = L_raw.flatten() / 100.0 # change into units R*/10ms
             
                M_raw = image['LMS_Image'][:,:,1]
                M_power[j,:] = sig.welch2d(M_raw[0:SIZE,0:SIZE])
                M_raw = M_raw.flatten() / 100.0

                S_raw = image['LMS_Image'][:,:,2]
                S_power[j,:] = sig.welch2d(S_raw[0:SIZE,0:SIZE])
                S_raw = S_raw.flatten() / 100.0
                
                L_count[j,:],bins = np.histogram(L_raw,BINS) 
                M_count[j,:],bins = np.histogram(M_raw,BINS)
                S_count[j,:],bins = np.histogram(S_raw,BINS)
                
                j += 1

        out = [{}]*4
        out[0] =   {
                     'name':'L',
                     'count':L_count,
                     'power':L_power
                     }
        out[1] =   {
                     'name':'M',
                     'count':M_count,
                     'power':M_power
                     }
        out[2] = {
                     'name':'S',
                     'count':S_count,
                     'power':S_power
                     }
        out[3] =     {
                    'bins': BINS
                    }             

    if Analysis_Option == 3:
    
        SIZE = sio.loadmat(files[0])
        SIZE = SIZE['LMS_Image'][:,:,0]
        SIZE = len(SIZE[:,0])
        L_power = np.zeros((np.sum(goodFiles), np.floor(SIZE / 2 + 1)))
        M_power = np.zeros((np.sum(goodFiles), np.floor(SIZE / 2 + 1)))
        S_power = np.zeros((np.sum(goodFiles), np.floor(SIZE / 2 + 1)))
        L_corr = np.zeros((np.sum(goodFiles),323))
        M_corr = np.zeros((np.sum(goodFiles),323))
        S_corr = np.zeros((np.sum(goodFiles),323))

        for i in range(0,stop):
            if goodFiles(i) == 1:
            
                image = sio.loadmat(files[i]) 

                L_raw = image['LMS_Image'][:,:,0]
                L_proc = sig.RemoveMean(L_raw)
                L_corr[j,:] = sig.pxl_corr(L_proc[100:422,100:422])
                L_power[j,:] = sig.welch2d(L_raw[0:SIZE,0:SIZE])
                L_raw = L_raw.flatten() / 100.0 # change into units R*/10ms
             
                M_raw = image['LMS_Image'][:,:,1]
                M_proc = sig.RemoveMean(M_raw)
                M_corr[j,:] = sig.pxl_corr(M_proc[100:422,100:422])
                M_power[j,:] = sig.welch2d(M_raw[0:SIZE,0:SIZE])
                M_raw = M_raw.flatten() / 100.0

                S_raw = image['LMS_Image'][:,:,2]
                S_proc = sig.RemoveMean(S_raw)
                S_corr[j,:] = sig.pxl_corr(S_proc[100:422,100:422])
                S_power[j,:] = sig.welch2d(S_raw[0:SIZE,0:SIZE])
                
                S_raw = S_raw.flatten() / 100.0
                
                L_count[j,:],bins = np.histogram(L_raw,BINS) 
                M_count[j,:],bins = np.histogram(M_raw,BINS)
                S_count[j,:],bins = np.histogram(S_raw,BINS)
                
                j += 1
        
        out = [{}]*4
        out[0] =   {
                     'name':'L',
                     'count':L_count,
                     'power':L_power,
                     'corr':L_corr
                     }
        out[1] =   {
                     'name':'M',
                     'count':M_count,
                     'power':M_power,
                     'corr':M_corr
                     }
        out[2] = {
                     'name':'S',
                     'count':S_count,
                     'power':S_power,
                     'corr':S_corr
                     }
        out[3] =     {
                    'bins': BINS
                    }             
                     
    return out         
    
    
### Information functions ###
#############################

def ArrayInfoFunc(XY_array,interp_method = 'poly',order = 2, res = 1000):
    """
    ArrayInfoFunc(P_ratio_XY,P_ratio_LS,Numer_XY,Numer_LS)
    
    Compute the information in an array of cones based on small array
    information rates computed by ConeInfoFunc
    
    :param XY_array: takes XY_array generated from ConeInfoFunc and generates
                    approximation of larger cone array using equation 3 from 
                    page 6 of garrigan et al.
    :param interp_method: 'poly' or 'linear'. Default = 'poly'.
    

    :returns: updated XY_array.
    
    .. note:: 
       Used after a call to ConeInfoFunc()
 
    """


    # First find the P_ratio_XY, i.e. Mutual/Total Information.
    XY_array['P_ratio_XY'] = XY_array['Mutual_XY'] / np.mean(XY_array['Info_XY'], axis = 0)

    # Second find the numerator to equation 3 on page 6 of Garrigan et al.
    XY_array['Numerator_XY'] = np.mean(XY_array['Info_XY'], axis = 0) - XY_array['Mutual_XY']



    if interp_method == 'linear': # linear interpolation to find P_XY:
        raise ValueError('linear interp not currently supported')
            ## add linear interoplation method?

            
    if interp_method == 'poly': # this appears to be something like what Garrigan et al. use.
        Pfit_XY = np.polyfit(np.arange(0,7),XY_array['P_ratio_XY'],order)
        P_XY = np.polyval(Pfit_XY,np.linspace(0,6,res))
        
        Nfit_XY = np.polyfit(np.arange(0,7),XY_array['Numerator_XY'],order)
        Num_XY = np.polyval(Nfit_XY,np.linspace(0,6,res))
    else:
        raise ValueError('Unexpected interpolation method')



    # Fourth find the Array Information using equation 3, page 6.
    Array_Info_XY = np.zeros((res))
    for i in range(0, res):
        Array_Info_XY[i] = Num_XY[i] / (1.0 + P_XY[i])


    XY_array['Array_Info_XY'] = Array_Info_XY / (6.0 * 10.0)
    return XY_array


def ConeInfoFunc(Dict,Cone_X, Cone_Y, Cone_Spacing = 1,Number_Cones = 6):
    """
    ConeInfoFunc(X_analysis,Y_analysis,Cone_Spacing, Number_Cones)
    
    Using the equations from Garrigan et al. for small arrays
    Computes: Info = Info_X + Info_Y - Info_XY  
    
    :param Dict: dict containing analysis files generated from SignalGen().
    :param Cone_X:
    :param Cone_Y:
    :param Number_Cones: number of cones in the array.
    :param Cones_Spacing: distance in pixels between cones in the array.
    
    :returns: outputs a dict with Information data and metadata.

    .. note::  
       * No options right now.  
       * Eventually want to derive Scale_Factor from newly collected data.
       * Called after SignalGen()
 
    """


    Num = np.arange(0, Number_Cones + 1) # number of X cones in array.
    
    X_analysis = Dict[Cone_X]
    Y_analysis = Dict[Cone_Y]

    X_count = X_analysis['count']
    Y_count = Y_analysis['count']
    
    # small 6 pixel arrays:
    I_XY = np.zeros(((X_count.shape[0]),7))
    Mutual_XY = np.zeros((7))


    Scale_Factor = [0.75, 0.84, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91] # from Garrigan et al.
     
    for i in range(0,7):
        Number_X = Num[i]

        ## Set up scaling params

        Number_Y = Number_Cones - Number_X
        

        # find distance with linear interpolation:
        if Number_X != 0:
            Distance_X = np.sqrt(Number_Cones * Cone_Spacing / Number_X)
            Delta_X = np.interp(Distance_X - 1.0, np.arange(0,len(Scale_Factor)), Scale_Factor)

        if Number_Y != 0:
            Distance_Y = np.sqrt(Number_Cones * Cone_Spacing / Number_Y)
            Delta_Y = np.interp(Distance_Y - 1.0, np.arange(0,len(Scale_Factor)), Scale_Factor)

            
        ## Information calculations:

        # compute the variance and information for each image serially, for each cone array type
        Info_X = np.zeros((X_count.shape[0]))
        Info_Y = np.zeros((Y_count.shape[0]))

        for j in range(0,X_count.shape[0]):

           # Information in each array:
            if Number_X != 0:
                signal_X = np.var(X_count[j,:])
                noise_X = np.mean(X_count[j,:])

                Info_X[j] = (0.5 * (np.log2(1.0 + signal_X / noise_X))) * Number_X**Delta_X
            elif Number_X == 0:
                Info_X[j] = 0
            

            if Number_Y != 0:
                signal_Y = np.var(Y_count[j,:])
                noise_Y = np.mean(Y_count[j,:])

                Info_Y[j] = (0.5 * (np.log2(1.0 + signal_Y / noise_Y))) * Number_Y**Delta_Y
            elif Number_Y == 0: 
                Info_Y[j] = 0

                
        # Setup for probability calculations
        if Number_X and Number_Y != 0:
            n = 16 
            Pbins = np.linspace(min(min(Info_X),min(Info_Y)),
                                max(max(Info_X),max(Info_Y)),n)
                                

            # Calculate joint distribution:
            JointProb,f,f = np.histogram2d(Info_X,Info_Y,Pbins) # count
            JointProb = (JointProb + 1.0) / np.sum(JointProb + 1.0) # prob distribution

            #Calculate marginal distribtions:
            Marg_X = np.sum(JointProb,axis=1)
            Marg_Y = np.sum(JointProb,axis=0)

            #Calculate mutual info between arrays:
            Info_XY = 0
            for k in range(0,len(Marg_X)):
                for l in range(0,len(Marg_Y)):
                    Info_XY += JointProb[k,l] * np.log2(Marg_X[k] * Marg_Y[l] / JointProb[k,l])
                
            #Invert sign
            Info_XY = -Info_XY

            #Calculate info in the array:
            Info = Info_X + Info_Y - Info_XY
        
        # Deal with zero cases:
        elif Number_X == 0:
            Info = Info_Y
            Info_XY = 0
        elif Number_Y == 0:
            Info = Info_X
            Info_XY = 0
        
        
        I_XY[:,i] = Info
        Mutual_XY[i] = Info_XY
        
    XY_array = {
                'Info_XY':I_XY,
                'Mutual_XY':Mutual_XY,
                'P_ratio_XY':[],
                'Numerator_XY':[],
                'Array_Info_XY':[],
                'Num':Num,
                'X_name':X_analysis['name'],
                'Y_name':Y_analysis['name']
                }
    return XY_array

    

def OptPercentX(XY_array,PrintOption = 1):
    """
    OptPercentX(XY_array,PrintOption)
    
    Find the optimal percentage of X in an XY cone array with respect to
    information.
    
    :param XY_array: structure containing data about XY cone array
    :param PrintOption: print result. Yes = 1, No = 0. 
    
    :returns: updated structure with OptPercentX
    

    .. note:: Called after ConeInfoFunc and ArrayInfoFunc
 
    """


    indXY = float(np.argmax(XY_array['Array_Info_XY']))
    XY_array['optPercentXY'] = indXY / len(XY_array['Array_Info_XY']) * 100.0

    if PrintOption == 1:
    
        print 'Optimal %',XY_array['X_name'],'in',XY_array['X_name'],',',XY_array['Y_name'], ':',str(XY_array['optPercentXY'])

    
    return XY_array



    

### Plotting funcs ###
######################

def SignalPlot(Dict,FONTSIZE = 22):
    """
 
    SignalPlots(L_analysis,M_analysis,S_analysis,BINS,Compute_Analyses)
    
    Plot statistical data about Signals: histograms, spatial frequency and
    average correlation.
    

    :param Dict: input to plot
    :param FONTSIZE: figure text fontsize
    
    : returns: Signal plot.

    
    """
    fig = {}
    
    if 'count' in Dict[0]:
    
        fig1 = plt.figure()
        
        font = {'weight' : 'norm', 'size'  : FONTSIZE}
        legend = {'frameon' : False}
        ticks = {'direction' : 'out', 'major.size' : 10 }
        #axes = {'limits'}
        plt.rc('font', **font)
        plt.rc('legend', **legend)
        plt.rc('xtick', **ticks)
        plt.rc('ytick', **ticks)
        
        ax = fig1.add_subplot(1, 1, 1)
        
        ax.plot(Dict[3]['bins'][1:], np.mean(Dict[0]['count'],0),'r',
                Dict[3]['bins'][1:], np.mean(Dict[1]['count'],0),'g',
                Dict[3]['bins'][1:], np.mean(Dict[2]['count'],0),'b', linewidth = 2.5)
                
        
        ax.legend(('L','M','S'))
        ax = pf.TufteAxis(ax,['left','bottom'])
        
        # force scientific notation:    
        #out = pf.SciNoteAxis(plt.gca(),['x','y'])
        
        #plt.gca().xaxis.major.formatter.set_powerlimits((-1, 0))
        #plt.axes.ticklabel_format(style='sci', axis = 'y') 
        #plt.gca().yaxis.major.formatter.set_powerlimits((-1, 0))
        #plt.axes.ticklabel_format(style='sci', axis = 'x') 
        
        plt.xlabel('R*/10ms')
        plt.ylabel('pixels')
        plt.tight_layout()
        plt.show(block = False)
        
        fig['fig1'] = fig1
        
    # Plot average spatial frequency
    if 'power' in Dict[0]:
    
        fig2 = plt.figure()
        font = {'weight' : 'norm', 'size'  : FONTSIZE}
        legend = {'frameon' : False}
        ticks = {'direction' : 'out' , 'major.size' : 10}
        plt.rc('font', **font)
        plt.rc('legend', **legend)
        plt.rc('xtick', **ticks)
        plt.rc('ytick', **ticks)
        ax = fig2.add_subplot(111)

        ax.semilogx(np.arange(1,Dict[0]['power'].shape[1] + 1)*2. / 46. / 2.,np.mean(Dict[0]['power'],0),
                    'r.', markersize = 10, markeredgewidth = 0)
        plt.hold(True)
        
        ax.semilogx(np.arange(1,Dict[0]['power'].shape[1] + 1)*2. / 46. / 2.,np.mean(Dict[1]['power'],0),
                    'g.', markersize = 10, markeredgewidth = 0)
        plt.hold(True)
        
        ax.semilogx(np.arange(1,Dict[0]['power'].shape[1] + 1)*2. /46. / 2.,np.mean(Dict[2]['power'],0),
                    'b.', markersize = 10, markeredgewidth = 0)
                    
        
        ax.legend(('L','M','S'))
        ax = pf.simpleaxis(ax)
        
        plt.xlabel('spatial frequency (cycles / deg)')
        plt.ylabel('spectral density (dB)')
        
        plt.tight_layout()
        plt.show(block = False)
        
        fig['fig2'] = fig2
        
    # Plot average correlation
    if 'corr' in Dict[0]:
    
        fig3 = plt.figure()
        font = {'weight' : 'norm', 'size'  : FONTSIZE}
        legend = {'frameon' : False}
        ticks = {'direction' : 'out', 'major.size' : 10 }
        plt.rc('font', **font)
        plt.rc('legend', **legend)
        plt.rc('xtick', **ticks)
        plt.rc('ytick', **ticks)
        ax = fig3.add_subplot(111)
        
        ax.semilogx(np.arange(1,len(Dict[0]['corr']) + 1) / 46., np.mean(Dict[0]['corr'],0),
                'r', linewidth = 2.5)
        plt.hold(True)
        
        ax.semilogx(np.arange(1,len(Dict[1]['corr']) + 1) / 46., np.mean(Dict[1]['corr'],0),
                'g', linewidth = 2.5)
        plt.hold(True)
        
        ax.semilogx(np.arange(1,len(Dict[2]['corr']) + 1) / 46., np.mean(Dict[2]['corr'],0),
                'b', linewidth = 2.5)
        
        ax = pf.TufteAxis(ax,['left','bottom'])
        ax.legend('L','M','S')
        #set(gca,'fontsize',25, 'linewidth',2, 'TickDir', 'out')
        ax.set_ylim([0, 1])
        ax.set_xlim([0, 7])
        ax.xlabel('seperation (deg)')
        ax.ylabel('mean correlation')
        plt.show(block = False)
        
        fig['fig3'] = fig3
        
    return fig

        
def InfoPlots(XY_array1,XY_array2, FONTSIZE = 22):
    """
    InfoPlots(XY_array1,XY_array2)
    
    Plot data Information data
    
    :param XY_array1: structure with Information from 2 cone array.
    :param XY_array2: structure with Information from 2 cone array.
    

    :returns: Plots

    .. note:: 
       Need to add options.
 
    """
    fig = {}
    ## plot Info results from arrays of different compositions:
    if 'Info_XY' in XY_array1:
        
        fig1 = plt.figure(1)
        font = {'weight' : 'norm', 'size'  : FONTSIZE}
        legend = {'frameon' : False}
        ticks = {'direction' : 'out' , 'major.size' : 10}
        plt.rc('font', **font)
        plt.rc('legend', **legend)
        plt.rc('xtick', **ticks)
        plt.rc('ytick', **ticks)
        ax = fig1.add_subplot(111)
        

        ax.plot(XY_array1['Num'] / max(XY_array1['Num']) * 100.0, np.mean(XY_array1['Info_XY'],0), 'k',
                XY_array2['Num'] / max(XY_array2['Num']) * 100.0, np.mean(XY_array2['Info_XY'],0), '--k', linewidth = 3)
        
        
        ax.legend(('L,M array','L,S array'), loc = 'lower center')
        ax = pf.TufteAxis(ax,['left','bottom'])
        plt.xlabel('% L')
        plt.ylabel('information (bits)')
        plt.tight_layout()
        plt.show(block=False)
        
        fig['fig1'] = fig1
        
    # plot mutual/info from %L
    if 'P_ratio_XY' in XY_array1:
    
        fig2 = plt.figure(2)
        font = {'weight' : 'norm', 'size'  : FONTSIZE}
        legend = {'frameon' : False}
        ticks = {'direction' : 'out' , 'major.size' : 10}
        plt.rc('font', **font)
        plt.rc('legend', **legend)
        plt.rc('xtick', **ticks)
        plt.rc('ytick', **ticks)
        ax = fig2.add_subplot(111)
        
        ax.plot(XY_array1['Num'] / max(XY_array1['Num']) * 100.0, XY_array1['P_ratio_XY'],'k',
                XY_array2['Num'] / max(XY_array2['Num']) * 100.0, XY_array2['P_ratio_XY'], '--k', linewidth = 3)

        ax.legend(('L,M array','L,S array'), loc = 'lower center')
        ax = pf.TufteAxis(ax,['left','bottom'])
        plt.xlabel('% L')
        plt.ylabel('mutual/total information')
        plt.tight_layout()
        plt.show(block=False)
        
        fig['fig2'] = fig2

    # plot optimal mosaic from %L
    if 'Array_Info_XY' in XY_array1:
    
        fig3 = plt.figure(3)
        font = {'weight' : 'norm', 'size'  : FONTSIZE}
        legend = {'frameon' : False}
        ticks = {'direction' : 'out' , 'major.size' : 10}
        plt.rc('font', **font)
        plt.rc('legend', **legend)
        plt.rc('xtick', **ticks)
        plt.rc('ytick', **ticks)
        ax = fig3.add_subplot(111)
        

        ax.plot(np.linspace(0,100,1000), XY_array1['Array_Info_XY'],'k',
                np.linspace(0,100,1000), XY_array2['Array_Info_XY'],'--k', linewidth = 3)

        ax.legend(('L,M array','L,S array'), loc = 'lower center')
        ax = pf.TufteAxis(ax,['left','bottom'])
        plt.xlabel('% L')
        plt.ylabel('information rate'); # (bits/(cone*10ms)
        plt.tight_layout()
        plt.show(block=False)
        
        fig['fig3'] = fig3
        
        return fig

        
def main():
    """
    
    Start here.
    
    Will get user input and then compute information.
    
    """
    
    # generate signal and plots:
    Signal = SignalGen()
    fig = SignalPlot(Signal)
    
    print ''
    
    for i,key in enumerate(fig):
    
        plt.figure(i + 1)
        
        save_option = raw_input("Save figure {0} (y / n)? ".format(i + 1))
            
        if save_option.lower() == 'y' or save_option.lower() == 'yes':
        
            plt.savefig('SignalPlot' + str(i + 1) + '.png')
            plt.close(i + 1)
            
        elif save_option.lower() == 'n' or save_option.lower() == 'no':
            
            plt.close(i + 1)
        
        else:
            
            print 'Sorry, answer not understood. This is a yes or no question. Please try again'
            
            save_option = raw_input("Save figure {0} (y / n)? ".format(i + 1))
            plt.savefig('SignalPlot' + str(i + 1) + '.png')
            plt.close(i + 1)
            
    print ''
    
    # compute information:
    LM_array = ConeInfoFunc(Signal,0,1)
    LS_array = ConeInfoFunc(Signal,1,2)
    
    LM_array = ArrayInfoFunc(LM_array)
    LS_array = ArrayInfoFunc(LS_array)
    
    LM_array = OptPercentX(LM_array)
    LS_array = OptPercentX(LS_array)
    
    print ''
    
    fig = InfoPlots(LM_array,LS_array)

    
    for i,key in enumerate(fig):
    
        plt.figure(i + 1)

        save_option = raw_input("Save figure {0} (y / n)? ".format(i + 1))    
        
        if save_option.lower() == 'y' or save_option.lower() == 'yes':
            
            plt.savefig('InfoPlot' + str(i + 2) + '.png')
            plt.close(i + 1)
            
        elif save_option.lower() == 'n' or save_option.lower() =='no':
            
            plt.close(i + 1)
        
        else:
            
            print 'Sorry, answer not understood. Please try again'
            save_option = raw_input("Save figure {0} (y / n)? ".format(i + 1))
            plt.savefig('InfoPlot' + str(i + 1) + '.png')
            plt.close(i + 1)
            
    print ''
    
    raw_input("Press ENTER to exit") 
 
        
### Statistics ###
##################

def KL_Diverg(P, Q, Dim=0, BIN_RESOLUTION=50, INDEPENDENCE_TEST=0,FUNCTION = 'variance',PRINT=1):
    """
    Computes the Kullback-Leibler Divergence between two probability distributions
    KL_Diverg = Sum( P * log2( P/ Q) )

    for more info,
    see: http://en.wikipedia.org/wiki/Kullback%E2%80%93Leibler_divergence
    and: http://www.snl.salk.edu/~shlens/kl.pdf 
    
    :param P: The observed distribution.
    :type P: np.array
    :param Q: The model distribution you are testing
    :type Q: np.array
    :param Dim: the dimension of data you are look at.
    :type Dim: int
    :param BIN_RESOLUTION: the bin size for probability calculation.
    :type BIN_RESOLUTION: int
    :param INDEPENDENCE_TEST: if 1 test independence of two variables.
    :type INDEPENDENCE_TEST: int
    :param FUNCTION: Decide which moment to use when first processing the data: 'mean' or 'variance'.
    :type FUNCTION: str
    :param PRINT: option of printing the results to the console, 0 = No, 1 = print
    :type PRINT: int
    
    :returns: KL-Divergence in bits (log2)
    :rtype: float

    .. warning:: 
       not fully tested.

    """
    if FUNCTION == 'mean':
        Func = np.mean
    elif FUNCTION == 'variance':
        Func = np.var

    if len(P.shape) > 1:
        P = P[:,Dim]
    elif len(P.shape) == 1:
        P = P[:]

    if len(Q.shape) > 1:
        Q = Q[:,Dim]
    elif len(P.shape) == 1:
        Q = Q[:]
    
    Ptemp = np.zeros(len(P)/10.)
    for i in range(0,int(len(P)/10.)):
        Ptemp[i] = Func(P[10.*i:(10.*i)+10.])

    Qtemp = np.zeros(len(Q)/10.)
    for i in range(0,int(len(Q)/10.)):
        Qtemp[i] = Func(Q[10.*i:(10.*i)+10.])
        
    
    BINS = np.linspace(min( min(Ptemp),min(Qtemp) ),max (max(Ptemp),max(Qtemp)),BIN_RESOLUTION)
                
    P,b = np.histogram(Ptemp,BINS)
    Q,b = np.histogram(Qtemp,BINS)


    # avoid divide by zeros
    Q = Q + 1.0
    P = P + 1.0

    # normalizing the P and Qs    
    Q = Q / np.sum(Q)
    P = P / np.sum(P)
        
    # decide whether computing standard KL-diverg or KL-diverge to assess independence of two variable.
    if INDEPENDENCE_TEST == 0:

        temp =  P*np.log2(P/Q);

        # decide if one or two dimensional case:
        if Q.shape[0] == 1:
            dist = np.sum(temp,axis=1)
        elif Q.shape[0] == P.shape[0]:
            dist = np.sum(temp)

    if INDEPENDENCE_TEST == 1: # testing independence

        JointDist = np.outer(Q,P)
        MargQ = np.sum(JointDist,axis=0)
        MargP = np.sum(JointDist,axis=1)
        MargQP = np.outer(MargQ,MargP)
        
        temp = JointDist*np.log2(JointDist/MargQP)
        dist = np.sum(temp)

    if PRINT == 1:
        print 'KL-diverg = ', dist

    return dist    
        
        
if __name__ == '__main__':
    
    main()
