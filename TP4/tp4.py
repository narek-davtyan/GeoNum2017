#------------------------------------------------------
#
#  TP4 : Subdivision curves
#  http://tiborstanko.sk/teaching/geo-num-2017/tp4.html
#  [03-Mar-2017]
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
#    DataPts  :  d x 2 matrix of control points
#
def ReadDatapoints( filename ) :
    datafile = open(filename,'r');
    # datapoints
    p, closed = np.fromstring(datafile.readline(),sep=' ',dtype=int)
    isClosed = closed==1
    DataPts = np.fromfile(datafile,count=2*p,sep=' ',dtype=float)
    DataPts = DataPts.reshape(-1,2)
    return DataPts


#-------------------------------------------------
# CHAIKIN()
# Perform one iteration of the corner cutting
# subdivision scheme for a closed polygon.
#
# Input
#    X0 :  n x 2 matrix, represents a polygon
#
# Output
#    X1 :  n x 2 matrix, subdivided polygon
#
def Chaikin( X0 ) :
    
    # number of points
    n = X0.shape[0]
    
    # upsample
    X1 = np.zeros([2*n,2])
    k=0
    for i in range(0,n) :
        X1[k,:] = 3.0/4.0*X0[(i%n),:] + 1.0/4.0*X0[((i+1)%n),:]
        k+=1
        X1[k,:] = 1.0/4.0*X0[(i%n),:] + 3.0/4.0*X0[((i+1)%n),:]
        k+=1
    return X1


#-------------------------------------------------
# CORNERCUTTING()
# Perform one iteration of the corner cutting
# subdivision scheme for a closed polygon.
#
# Input
#    X0 :  n x 2 matrix, represents a polygon
#   a,b :  corner cutting parameters
#
# Output
#    X1 :  n x 2 matrix, subdivided polygon
#
def CornerCutting( X0, a, b ) :
    
    # number of points
    n = X0.shape[0]
    
    # upsample
    X1 = np.zeros([2*n,2])        
    k=0
    for i in range(0,n) :
        X1[k,:] = (1-a)*X0[(i%n),:] + a*X0[((i+1)%n),:]
        k+=1
        X1[k,:] = (1-b)*X0[(i%n),:] + b*X0[((i+1)%n),:]
        k+=1
    
    return X1


#-------------------------------------------------
# FOUTPOINT()
# Perform one iteration of the four-point
# subdivision scheme for a closed polygon.
#
# Input
#    X0 :  n x 2 matrix, represents a polygon
#     w :  tension parameter (generalized four-point only)
#
# Output
#    X1 :  n x 2 matrix, subdivided polygon
#
def FourPoint( X0, w ) :
    
    # number of points
    n = X0.shape[0]
    
    # upsample
    X1 = np.zeros([2*n,2])
    k=0
    for i in range(0,n) :
        X1[k,:] = X0[(i%n),:]
        k+=1
        X1[k,:] = -w*X0[((i-1)%n),:] + (1.0/2.0+w)*X0[(i%n),:] + \
                    (1.0/2.0+w)*X0[((i+1)%n),:] - w*X0[((i+2)%n),:]
        k+=1    
    return X1


#-------------------------------------------------
if __name__ == "__main__":
    
    # helper function for getting scheme name from ID
    def fullname(x):
        return {
            "CH" : "Chaikin",
            "CC" : "Corner cutting",
            "FP" : "Four-point",
        }.get(x,"Invalid")

    ###############################
    ## arg 1 : data name 
    ###############################
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
    else :
        dataname = "simple" # [simple,infinity]
    ###############################
    ## arg 2 : subdivision scheme
    ###############################
    if len(sys.argv) > 2 :
        scheme = sys.argv[2]
    else :
        scheme = "CH"
    if fullname(scheme) == "Invalid" :
        print " error :  invalid scheme "+scheme
        print "          should be CH, CC or FP"
        sys.exit(0)
    ###############################
    ## arg 3 : depth of subdivision
    ###############################
    if len(sys.argv) > 3 :
        depth = int(sys.argv[3])
    else :
        depth = 3
    
    # output : scheme name and subdivision depth
    print " "+fullname(scheme)
    print " depth = " + str(depth)
    
    # filename
    filename = DATADIR + dataname + ".data"
    
    # check if valid datafile
    if not os.path.isfile(filename) :
        print " error :  invalid dataname '" + dataname + "'"
        print " usage :  python tp4.py  [data=simple,infinity,bone,bunny]  [scheme=CH,CC,FP]  [depth=3]"
        
    else :

        ## (a,b) = (0.1,0.6) (b = a + 1/2)
        ## we obtain a C1 curve
        ## (a,b) = (0.1,0.5) (b != a + 1/2)
        ## we obtain a C0 curve 
        ## -- Corner cutting
        a = 0.1
        b = 0.5
        ## -- Generalized four-point
        ## W = (sqrt(5)-1)/8  => we obtain a fractal curve
        ## W(k) = W/k for k = 1...infinity => the curve converges to the polygon 
        w = (np.sqrt(5)+1)/8.0 + (np.sqrt(5)-1)/16.0
        
        # read datapoints
        DataPts = ReadDatapoints(filename)    
        
        # init subdivided curve
        SubPts = DataPts
        
        # iterative refinement
        for iteration in range(depth) :
            
            # Chaikin
            if scheme == "CH" :
                print(SubPts)
                SubPts = Chaikin(SubPts)
                print(SubPts)
            # Corner cutting
            elif scheme == "CC" :
                SubPts = CornerCutting(SubPts,a,b)
            
            # Four-point
            else :
                SubPts = FourPoint(SubPts,w)
           
        # set axes with equal proportions
        plt.axis('equal')
        
        # clear plot
        plt.cla()
        
        # plot coarse polygon
        plt.fill( DataPts[:,0], DataPts[:,1], edgecolor=.33*np.ones(3), linewidth=1, linestyle='--', fill=False)
        
        # plot subdivided polygon
        plt.fill( SubPts[:,0], SubPts[:,1], edgecolor='b', linestyle='-', linewidth=2, fill=False ) 
            
        # titles
        plt.gcf().canvas.set_window_title('TP4 Subdivision curves')
        plt.title(dataname+': scheme='+fullname(scheme)+', depth='+str(depth))

        ##
        ## TODO : Uncomment if you want to save the render as png image in data/
        ##
        #plt.savefig( DATADIR + dataname + ".png" )
        
        # render
        plt.show()
        