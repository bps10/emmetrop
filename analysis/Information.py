from __future__ import division
import numpy as np



def LargeArrayInfoFunc(XY_array,interp_method = 'poly',order = 2, res = 1000):
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
    XY_array['P_ratio_XY'] = XY_array['Mutual_XY'] / np.mean(XY_array['Info_XY'],
                                                             axis = 0)

    # Second find the numerator to equation 3 on page 6 of Garrigan et al.
    XY_array['Numerator_XY'] = (np.mean(XY_array['Info_XY'], axis = 0) - 
                                XY_array['Mutual_XY'])



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


def getScaleFactors():
    """Get the scale factors used in ConeInfoFunc.
    
    Currently estimated from Garrigan et al.  Eventually want to be able to
    rederive this from new image data.
    
    No input parameters
    
    :returns: approximated scale factors from Garrigan et al.
    
    .. todo::
       * Eventually want to derive Scale_Factor from newly collected data.
       
    """
    # from Garrigan et al.
    scale_factor = [0.75, 0.84, 0.86, 0.87, 0.88, 0.89, 0.90, 0.91]

    return scale_factor
    
def findConeSpacing(Number_X, Cone_Spacing, Scale_Factor, Number_Cones):
    """Find the distance between cones in a two cone type (X,Y) mosaic
    
    :param Number_X: number of X cones in array.
    :param Cone_Spacing: spacing between all cones.
    :param Scale_Factor: of the information.
    :param Number_Cones: number of cones total in the array.
    
    :returns: * Delta_X - correlation factor
              * Delta_Y - correlation factor
    """
    
    ## Set up scaling params

    Number_Y = Number_Cones - Number_X
    

    # find distance with linear interpolation:
    if Number_X != 0:
        Distance_X = np.sqrt(Number_Cones * Cone_Spacing / Number_X)
        Delta_X = np.interp(Distance_X - 1.0, 
                            np.arange(0,len(Scale_Factor)), Scale_Factor)

    if Number_Y != 0:
        Distance_Y = np.sqrt(Number_Cones * Cone_Spacing / Number_Y)
        Delta_Y = np.interp(Distance_Y - 1.0, 
                            np.arange(0,len(Scale_Factor)), Scale_Factor)
                                
    return Delta_X, Delta_Y
    

def findConeInfo(ConeNumber, ConeHist, Delta):
    """Find the information in a single cone type.
    
    :param ConeNumber: number of cone type
    :param ConeHist: 1D array of histogram data.
    :param Delta: correlation factor, based on cone spacing, derived from data.
    
    :returns: Info
    
    """
    
    if ConeNumber == 0:
        Info = 0
        
    elif ConeNumber != 0:
        signal_X = np.var(ConeHist)
        noise_X = np.mean(ConeHist)
        Info = (0.5 * (np.log2(1.0 + signal_X / noise_X))) * ConeNumber**Delta 
        
    return Info



def findMutualInfo(Number_X, Number_Y, Info_X, Info_Y):
    """Find the mutual information between two arrays
    
    Computes information according to: Info = Info_X + Info_Y - Mutual_Info    
    
    :param Number_X: number of X cones in an array.
    :param Number_Y: number of Y cones in an array.
    :param Info_X: information in X, computed from findConeInfo()
    :param Info_Y: information in Y, computed from findConeInfo()
    
    :returns: * Info - after accounting for mutual info
              * Mutual_Info - between the arrays
    """
    
    # Setup for probability calculations
    if Number_X and Number_Y != 0:
        n = 16 
        Pbins = np.linspace(min(min(Info_X),min(Info_Y)),
                            max(max(Info_X),max(Info_Y)),n)
                            

        # Calculate joint distribution:
            # count
        JointProb,f,f = np.histogram2d(Info_X,Info_Y,Pbins) 
            # prob distribution
        JointProb = (JointProb + 1.0) / np.sum(JointProb + 1.0) 

        #Calculate marginal distribtions:
        Marg_X = np.sum(JointProb,axis=1)
        Marg_Y = np.sum(JointProb,axis=0)

        #Calculate mutual info between arrays:
        Mutual_Info = 0
        for k in range(0,len(Marg_X)):
            for l in range(0,len(Marg_Y)):
                Mutual_Info += JointProb[k,l] * np.log2(Marg_X[k] * Marg_Y[l] /
                                                    JointProb[k,l])
            
        #Invert sign
        Mutual_Info = - Mutual_Info

        #Calculate info in the array:
        Info = Info_X + Info_Y - Mutual_Info
    
    # Deal with zero cases:
    elif Number_X == 0:
        Info = Info_Y
        Mutual_Info = 0
    elif Number_Y == 0:
        Info = Info_X
        Mutual_Info = 0

    return Info, Mutual_Info     
    

def SingleConeEntropyFunc(Cone_Hist, ConeNumber = 6):
    """Compute the entropy (bits) in a single cone
    
    :param Cone_Hist:
    :type Cone_Hist: np.array
    :param ConeNumber: a list or int indicating the number(s) of cones in
                        simulated mosaic of cones.
    :returns: an array of entropy calculations in shape of ConeNumber.
    
    For a mixed mosaic of cones, use SmallArrayInfoFunc()
    
    """
    
    Scale_Factor = getScaleFactors()
    Entropy = np.zeros((len(ConeNumber)))
    for conenum in ConeNumber:
        
        index = conenum - 1 # subtract 1 to account for 0...n-1 indexing.
        #for i in range(0,Cone_Hist.shape[0]):
        
        Entropy[index] = findConeInfo(conenum, Cone_Hist, 
                                   Scale_Factor[index])
        
    return Entropy    

           
def SmallArrayInfoFunc(Cone_X_Hist, Cone_X_type, 
                       Cone_Y_Hist = None, Cone_Y_type = None, 
                       Cone_Spacing = 1, Number_Cones = 6):
    """Using the equations from Garrigan et al. for small arrays
    Computes: Info = Info_X + Info_Y - Info_XY  
    
    :param Cone_X_Hist: histogram data for group of LMS images
    :param Cone_X_type: cone type (LMS)
    :param Cone_Y_Hist: histogram data for group of LMS images
    :param Cone_Y_type: cone type (LMS)
    :param Cone_Spacing: distance in pixels between cones in the array.
    :param Number_Cones: number of cones in an array.
    
    :returns: outputs a dict with Information data and metadata.
 
    """
    Scale_Factor = getScaleFactors()
           
    Num = np.arange(0, Number_Cones + 1) # number of X cones in array.

    # small 6 pixel arrays:
    Info_XY = np.zeros(((Cone_X_Hist.shape[0]),7))
    Mutual_XY = np.zeros((7))

        
    for i in Num:
        
        Number_X = Num[i]
        Number_Y = Number_Cones - Number_X
        
        # Find cone spacing:
        Delta_X, Delta_Y = findConeSpacing(Number_X, Cone_Spacing, 
                                           Scale_Factor, Number_Cones)    
        ## Information calculations:
        # compute the variance and information for each image serially, 
        # for each cone array type
        Info_X = np.zeros((Cone_X_Hist.shape[0]))
        Info_Y = np.zeros((Cone_Y_Hist.shape[0]))

        for j in range(0,Cone_X_Hist.shape[0]):

           # Information in each array:

            Info_X[j] = findConeInfo(Number_X, Cone_X_Hist[j,:], Delta_X)
            Info_Y[j] = findConeInfo(Number_Y, Cone_Y_Hist[j,:], Delta_Y)
            

        Info, Mutual_Info = findMutualInfo(Number_X, Number_Y, 
                                           Info_X, Info_Y)
        
        Info_XY[:,i] = Info
        Mutual_XY[i] = Mutual_Info
        
    XY_array = {
                'Info_XY': Info_XY,
                'Mutual_XY': Mutual_XY,
                'Num': Num,
                'X_name': Cone_X_type,
                'Y_name': Cone_Y_type
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
    
        print 'Optimal %',XY_array['X_name'],'in',XY_array['X_name'],',', \
        XY_array['Y_name'], ':',str(XY_array['optPercentXY'])

    
    return XY_array

 
        
### Statistics ###
##################

def KL_Diverg(P, Q, Dim=0, BIN_RESOLUTION=50, INDEPENDENCE_TEST=0,
              FUNCTION = 'variance',PRINT=1):
    """
    Computes the Kullback-Leibler Divergence between two probability
    distributions

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
    :param FUNCTION: Decide which moment to use when first processing the \
                     data: 'mean' or 'variance'.
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
        
    
    BINS = np.linspace(min( min(Ptemp),min(Qtemp) ),max (max(Ptemp),
                           max(Qtemp)),BIN_RESOLUTION)
                
    P,b = np.histogram(Ptemp,BINS)
    Q,b = np.histogram(Qtemp,BINS)


    # avoid divide by zeros
    Q = Q + 1.0
    P = P + 1.0

    # normalizing the P and Qs    
    Q = Q / np.sum(Q)
    P = P / np.sum(P)
        
    # decide whether computing standard KL-diverg or KL-diverge to assess 
    # independence of two variable.
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

def main():
    
    raise('Sorry a direct call to information is not supported. See doc')        
        
if __name__ == '__main__':
    
    main()
