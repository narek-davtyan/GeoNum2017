#------------------------------------------------------
#
#  TP2 : Bezier splines, C^k smoothness
#  http://tiborstanko.sk/teaching/geo-num-2017/tp2.html
#  [10-Feb-2017]
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
# READPOLYGON()
# Read datapoints from a file.
#
# Input
#    filename :  file to be read
#
# Output
#    DataPts  :  (n+1) x 2 matrix of datapoints
#
def ReadData( filename ) :
    datafile = open(filename,'r');
    l = datafile.readline()
    degree = np.fromstring(l,sep=' ',dtype=int)
    DataPts = np.fromfile(datafile,count=2*(degree+1),sep=' ',dtype=float)
    DataPts = DataPts.reshape(-1,2)
    return DataPts


#-------------------------------------------------
# DECASTELJAU( ... )
# Perform the De Casteljau algorithm.
#
# Input
#    BezierPts :  (degree+1) x 2 matrix of Bezier control points
#    k         :  upper index of the computed point (depth of the algorithm)
#    i         :  lower index of the computed point
#    t         :  curve parameter in [0.0,1.0]
#
# Output
#    point b_i^k from the De Casteljau algorithm.
#
def DeCasteljau( BezierPts, k, i, t ) :
    if k==0 :
        return BezierPts[i,:]
    else :
        return (1-t) * DeCasteljau( BezierPts, k-1, i, t )  +  t * DeCasteljau( BezierPts, k-1, i+1, t )


#-------------------------------------------------
# BEZIERCURVE( ... )
# Compute points on the Bezier curve.
#
# Input
#    BezierPts :  (degree+1) x 2 matrix of Bezier control points
#    N         :  number of curve samples
#    
# Output
#    CurvePts  :  N x 2 matrix of curvepoints
#
def BezierCurve( BezierPts, N ) :
    degree = BezierPts.shape[0]-1
    CurvePts = np.empty([N,2])
    i=0
    for t in np.linspace(0.0, 1.0, num=N) :
        CurvePts[i,:] = DeCasteljau( BezierPts, degree, 0, t )
        i+=1
    return CurvePts


#-------------------------------------------------
# COMPUTESPLINEC1( ... )
# Compute Bezier control points for a C1 quadratic spline interpolating the given data.
#
# Input
#    DataPts   : (n+1) x 2 matrix of datapoints P_0,...,P_n
#    
# Output
#    BezierPts : (2n+1) x 2 matrix of Bezier control points 
#
def ComputeSplineC1( DataPts ) :
    
    # input: DataPts = [P_0; P_1; ... ; P_n]
    n = DataPts.shape[0]-1
    
    # BezierPts matrix will contain (2n+1) rows :
    #   n+1 rows for input data and n rows for inner Bezier control pts.
    BezierPts = np.zeros([2*n+1,2])

    # fill BezierPoints with DataPoints 
    j = 0
    for i in range(0,n+1) :
        BezierPts[j] = DataPts[i]
        j += 2

    # set intermediate point for the first polygon as 
    # half of the distance of two points
    BezierPts[1] = 0.5*(BezierPts[0]+BezierPts[2])

    # we can vary the value of the first intermediate point
    # BezierPts[1] = 0.3*(BezierPts[0]+BezierPts[2]) - 0.1*BezierPts[4]

    # calculate the rest of intermediate points b_1 for other segments
    for i in range(1,2*n-1,2) :
        BezierPts[i+2] = 2*BezierPts[i+1] - BezierPts[i]
    
    return BezierPts


#-------------------------------------------------
# COMPUTESPLINEC2( ... )
# Compute Bezier control points for a C2 cubic spline interpolating the given data.
#
# Input
#    DataPts   : (n+1) x 2 matrix of datapoints P_0,...,P_n
#    
# Output
#              : (3n+1) x 2 matrix of Bezier control points 
#
def ComputeSplineC2( DataPts ) :

    # input: DataPts = [P_0; P_1; ... ; P_n]
    n = DataPts.shape[0]-1
    
    # matrix of the system
    M = np.zeros([3*n+1,3*n+1])
    
    # right side
    R = np.zeros([3*n+1,2])
    
    ##
    ## TODO : Fill the matrix of the system M.
    ##          C0  :  n+1 rows
    ##          C1  :  n-1 rows
    ##          C2  :  n-1 rows
    ##          bd  :    2 rows (natural spline)
    ##
    ## with 4 datapoints and 10 Bezier points A ... J,
    ## the matrix has the following form:
    ## M =                                               : corresp. equation   [type]
    ##  |  1   0   0   0   0   0   0   0   0   0  |      : A = A                [C0]
    ##  |  0   0   0   1   0   0   0   0   0   0  |      : D = D                [C0]
    ##  |  0   0   0   0   0   0   1   0   0   0  |      : G = G                [C0]
    ##  |  0   0   0   0   0   0   0   0   0   1  |      : J = J                [C0]
    ##  |  0   0   1  -2   1   0   0   0   0   0  |      : C - 2D + E = 0       [C1]
    ##  |  0   0   0   0   0   1  -2   1   0   0  |      : F - 2G + H = 0       [C1]
    ##  |  0   1  -2   0   2  -1   0   0   0   0  |      : B - 2C + 2E - F = 0  [C2]
    ##  |  0   0   0   0   1  -2   0   2  -1   0  |      : E - 2F + 2H - I = 0  [C2]
    ##  |  1  -2   1   0   0   0   0   0   0   0  |      : A - 2B + C = 0       [natural]
    ##  |  0   0   0   0   0   0   0   1  -2   1  |      : H - 2I + J = 0       [natural]
    ##

    # C0 : indicate with 1s the location of DataPoints to ensure continuity C0
    j = 0
    for i in range(0,n+1) :
        M[i][j] = 1
        j += 3
    
    # C1 : indicate coefficients to ensure continuity C1
    j = 2
    for i in range(n+1,2*n) :
        M[i][j+0] = 1
        M[i][j+1] = -2
        M[i][j+2] = 1
        j += 3
    
    # C2 : indicate coefficients to ensure continuity C2
    j = 1
    for i in range(2*n,3*n-1) :
        M[i][j+0] = 1
        M[i][j+1] = -2
        M[i][j+3] = 2
        M[i][j+4] = -1
        j += 3
    
    # bd : indicate coefficients to ensure natural spline
    j = 0
    for i in range(3*n-1,3*n+1) :
        M[i][j+0] = 1
        M[i][j+1] = -2
        M[i][j+2] = 1
        j = 3*n-2

    # put DataPts to the first (n+1) rows of R
    R[0:n+1] = DataPts[:]

    # return the solution
    return np.linalg.solve(M, R)
    
#-------------------------------------------------
if __name__ == "__main__":
    
    # arg 1 : data name 
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
    else :
        # spiral example can't work properly with C1,
        # due to the fact that the last segments are too small
        # and the intermidiate points of the last polygons
        # are too distorted to be C1-smooth, so the curve 
        # changes direction at the end
        dataname = "infinity" # simple, infinity, spiral, semi, tuple

    # arg 2 : sampling density
    if len(sys.argv) > 2 :
        density = int(sys.argv[2])
    else :
        density = 10

    # arg 3 : C2 continuity
    if len(sys.argv) > 3 :
        c2 = True
    else :
        c2 = False
    
    # uncomment this for manual C2 mode activation
    # c2 = True
    
    # filename
    filename = DATADIR + dataname + ".bcv"
    
    # check if valid datafile
    if not os.path.isfile(filename) :
        print "error:  invalid dataname '" + dataname + "'"
        print "usage:  python tp2.py  [simple,infinity,semi,spiral,tuple]  [sampling_density]  [c2]"
        
    else :    
        # read points to be interpolated
        DataPts = ReadData(filename)
        n = DataPts.shape[0]-1

        # compute Bezier points
        if c2 :
            BezierPts = ComputeSplineC2( DataPts )
            cstr='C2'
            deg=3
        else :
            BezierPts = ComputeSplineC1( DataPts )
            cstr='C1'
            deg=2

        # for each segment : compute and plot
        for i in range(0,n*deg,deg) :

            # build each control polygon with appropriate dimensions
            iBezierPts = np.zeros([deg+1,2])
            iBezierPts[:] = BezierPts[i:i+deg+1]
            
            # compute ant plot each segment
            CurvePts = BezierCurve( iBezierPts, density )
            plt.plot( CurvePts[:,0], CurvePts[:,1], '-', linewidth=3 )


        # plot the datapoints
        plt.plot( DataPts[:,0], DataPts[:,1], 'k.', markersize=10 )
        
        # plot the control polygon
        plt.plot( BezierPts[:,0], BezierPts[:,1], 'k.--', linewidth=1 )
        
        # set axes with equal proportions
        plt.axis('equal')
        
        # titles
        plt.gcf().canvas.set_window_title('TP2 Bezier splines')
        plt.title(cstr+' '+dataname+', '+str(density)+" pts/segment")
        
        # save the render as png image in the data/ dir
        #plt.savefig( DATADIR + dataname + ".png" )
        
        # render
        plt.show()