#------------------------------------------------------
#
#  TP5 : Lane-Riesenfeld algorithm
#  http://tiborstanko.sk/teaching/geo-num-2017/tp5.html
#  [10-Mar-2017]
#
#------------------------------------------------------
#
#  This file is a part of the course:
#    Geometrie numerique (spring 2017)
#    https://github.com/GeoNumTP/GeoNum2017
#    M1 Informatique
#    UFR IM2AG
#
#  Course lecturer:
#    Georges-Pierre.Bonneau at inria.fr
#
#  Practical part:
#    Tibor.Stanko at inria.fr
#
#------------------------------------------------------

import sys, os
import matplotlib.pyplot as plt
import numpy as np

TP = os.path.dirname(os.path.realpath(__file__)) + "/"
DATADIR = filename = TP+"data/"


#-------------------------------------------------
# READDATAPOINTS()
# Read datapoints from a file.
#
# Input
#    filename :  file to be read
#
# Output
#    DataPts  :  d x 2 matrix, represents a closed polygon
#
def ReadDatapoints( filename ) :
    datafile = open(filename,'r');
    p = np.fromstring(datafile.readline(),sep=' ',dtype=int)
    DataPts = np.fromfile(datafile,count=2*p,sep=' ',dtype=float)
    DataPts = DataPts.reshape(-1,2)
    return DataPts


#-------------------------------------------------
# LANERIESENFELD()
# Perform one iteration of the Lane-Riesenfeld algorithm.
#
# Input
#    X0       :  n x 2 matrix, initial polygon
#    degree   :  degree of subdivision
#
# Output
#    X1       :  2n x 2 matrix, subdivided polygon
#
#
# reference paper:
#   Lane, Riesenfeld
#   A Theoretical Development for the Computer Generation
#    and Display of Piecewise Polynomial Surfaces
#   1980. IEEE Transactions on Pattern Analysis and Machine Intelligence
#   https://doi.org/10.1109/TPAMI.1980.4766968
#
#
def LaneRiesenfeld(X0,degree) :

    # number of points
    n = X0.shape[0]
    
    # upsample
    X1 = np.zeros([2*n,2])
    k=0

    # refining
    for i in range(0,n) :
        X1[k,:] = X0[i,:]
        k+=1
        X1[k,:] = X0[i,:]
        k+=1
    
    # smoothing
    for d in range(0,degree) :
        tmp = np.zeros([2*n,2])
        for i in range(0,k) :
            tmp[i,:] = 0.5*(X1[(i%k),:] + X1[((i+1)%k),:])
        
        # replacement
        X1 = tmp

    return X1

#-------------------------------------------------
# LANERIESENFELD2()
# Perform one iteration of the Lane-Riesenfeld algorithm for degree2.
#
# Input
#    X0       :  n x 2 matrix, initial polygon
#
# Output
#    X1       :  2n x 2 matrix, subdivided polygon
#
def LaneRiesenfeld2(X0) :

    # number of points
    n = X0.shape[0]
    
    # upsample + refining
    X1 = np.zeros([2*n,2])
    k=0
    for i in range(0,n) :
        X1[k,:] = 3.0/4.0*X0[(i%n),:] + 1.0/4.0*X0[((i+1)%n),:]
        k+=1
        X1[k,:] = 1.0/4.0*X0[(i%n),:] + 3.0/4.0*X0[((i+1)%n),:]
        k+=1

    return X1

#-------------------------------------------------
# FOURPOINT()
# Perform one iteration of the 4-point Lane-Riesenfeld subdivision.
# 
# Input
#    X0       :  n x 2 matrix, initial polygon
#    degree   :  degree of subdivision
#
# Output
#    X1       :  2n x 2 matrix, subdivided polygon
#
#
# reference paper:
#   Cashman et al.
#   Generalized Lane-Riesenfeld algorithms
#   2013. Computer Aided Geometric Design
#   http://dx.doi.org/10.1016/j.cagd.2013.02.001
#
#
def FourPoint(X0,degree) :

    # number of points
    n = X0.shape[0]
    #print("degree p = "+str(n))
    # upsample
    X1 = np.zeros([2*n,2])
    k=0

    # refining
    for i in range(0,n) :
        X1[k,:] = X0[i,:]
        k+=1
        X1[k,:] = 1.0/16.0*(-X0[(i-1)%n,:] + 9*X0[(i%n),:] \
            + 9*X0[((i+1)%n),:] - X0[((i+2)%n),:])
        k+=1
    
    # smoothing
    for d in range(0,degree) :
        tmp = np.zeros([2*n,2])
        for i in range(0,k) :
            tmp[i,:] = 1.0/16.0*(-X1[(i-1)%k,:] + 9*X1[i,:] \
                + 9*X1[((i+1)%k),:] - X1[((i+2)%k),:])
    
        # replacement
        X1 = tmp

    return X1


#-------------------------------------------------
# SIXPOINT()
# Perform one iteration of the 6-point Lane-Riesenfeld subdivision.
# 
# Input
#    X0       :  n x 2 matrix, initial polygon
#    degree   :  degree of subdivision
#
# Output
#    X1       :  2n x 2 matrix, subdivided polygon
#
#
# reference paper:
#   Ashraf et al.
#   A Six-Point Variant on the Lane-Riesenfeld Algorithm
#   2014. Journal of Applied Mathematics
#   http://dx.doi.org/10.1155/2014/628285
#
#
def SixPoint(X0,degree) :
    
    # number of points
    n = X0.shape[0]
    
    # upsample
    X1 = np.zeros([2*n,2])
    k=0

    # refining
    for i in range(0,n) :
        X1[k,:] = X0[i,:]
        k+=1
        X1[k,:] = 1.0/256.0*(3*X0[(i-2)%n,:] -25*X0[(i-1)%n,:] + 150*X0[(i%n),:] \
            + 150*X0[((i+1)%n),:] - 25*X0[((i+2)%n),:] + 3*X0[((i+3)%n),:])
        k+=1
    
    # smoothing
    for d in range(0,degree) :
        tmp = np.zeros([2*n,2])
        for i in range(0,k) :
            tmp[i,:] = 1.0/256.0*(3*X1[(i-2)%k,:] -25*X1[(i-1)%k,:] + 150*X1[(i%k),:] \
                + 150*X1[((i+1)%k),:] - 25*X1[((i+2)%k),:] + 3*X1[((i+3)%k),:])    
    
        # replacement
        X1 = tmp

    return X1


#-------------------------------------------------
if __name__ == "__main__":
    
    #---------------------------------------------
    # arg 1 : data name 
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
    else :
        dataname = "hepta" #bone,infinity,sumsign
        
    #---------------------------------------------    
    # arg 2 : scheme name
    if len(sys.argv) > 2 :
        scheme = sys.argv[2]
    else :
        scheme = "LR"
        
    # check
    if scheme == "LR" :
        schemeName = "Lane-Riesenfeld"
    elif scheme == "LR2" :
        schemeName = "Lane-Riesenfeld2"
    elif scheme == "FP" :
        schemeName = "4-point"
    elif scheme == "SP" :
        schemeName = "6-point"
    else :            
        print " error :  invalid scheme "+scheme
        print "          should be LR or LR2 or FP or SP"
        sys.exit(0)
        
    #---------------------------------------------
    # arg 3 : degree of the curve
    if len(sys.argv) > 3 :
        degree = int(sys.argv[3])
    else :
        degree = 2
        
    #---------------------------------------------
    # arg 4 : number of subdivisions
    if len(sys.argv) > 4 :
        subdivisions = int(sys.argv[4])
    else :
        subdivisions = 5
        
    # filename
    filename = DATADIR + dataname + ".data"
    
    # check if valid datafile
    if not os.path.isfile(filename) :
        print " error :  invalid dataname '" + dataname + "'"
        print " usage :  python tp5.py  [data=hepta; bone,infinity,sumsign]  [scheme=LR; FP,SP]  [degree=3]  [subdivisions=5]"
        
    else :
        
        # read data
        P = ReadDatapoints(filename)    
        
        # init subdivided polygon
        X = P
        #print(X.shape[0])
        # iterative subdivision
        for i in range(subdivisions) :
            
            # Lane-Riesenfeld
            if scheme == "LR" :
                X = LaneRiesenfeld(X,degree)

            # Lane-Riesenfeld 2
            if scheme == "LR2" :
                X = LaneRiesenfeld2(X)

            # 4-point
            elif scheme == "FP" :
                X = FourPoint(X,degree)
                
            # 6-point
            else :
                X = SixPoint(X,degree)
        
        print(X)
        # set axes with equal proportions
        plt.axis('equal')
        
        # plot the data
        plt.cla()
        
        # input polygon
        plt.fill( P[:,0], P[:,1], edgecolor=.33*np.ones(3), linewidth=1, linestyle='--', fill=False)
    
        # refined polygon
        plt.fill( X[:,0], X[:,1], edgecolor='b', linestyle='-', linewidth=2, fill=False )
            
        # titles
        # plot
        ptitle  = dataname+" : "
        ptitle += schemeName +", "
        ptitle += "deg="+str(degree)+", "
        ptitle += "sub="+str(subdivisions)
        plt.title(ptitle)
        # figure
        plt.gcf().canvas.set_window_title('TP5 Lane-Riesenfeld')

        ##
        ## TODO : Uncomment if you want to save the render as png image in data/ folder
        ##
        #plt.savefig( DATADIR+dataname+".png" )
        
        # render
        plt.show()        
        ## 2. When varying the degree with the three subdivision schemes:
        ##    LR: we observe that the more we increase the degree the more
        ##        the curve gets further from the control polygon, the result curve is Ck-1
        ##    FP: we observe no changes with the different datasets if we only slightly
        ##        change the degree, to observe the important change of the curve
        ##        on has to heavily change the degree parameter (e.g. with infinity - 30),
        ##        the result curve is C1
        ##    SP: we observe no changes with the different datasets if we only slightly
        ##        change the degree, to observe the important change of the curve
        ##        on has to even more heavily change the degree parameter, than with
        ##        the FP scheme (e.g. with infinity - 60),
        ##        the result curve is C2
        ## 3. The six-point scheme gives the best result in terms of smoothness
        ##    right away, that is with small degree (C2 curve), while the LR algorithm 
        ##    can get to the more substantial curve much faster (one can see the important 
        ##    change every time the degree parameter changes by one), the FP scheme
        ##    is in the middle, although its properties are much closer to the SP scheme.
        