#------------------------------------------------------
#
#  TP1 : Bezier curves, De Casteljau's algorithm
#  http://tiborstanko.sk/teaching/geo-num-2017/tp1.html
#  [03-Feb-2017]
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

TP = os.path.dirname(os.path.realpath(__file__)) + "\\"
DATADIR = filename = TP+"data\\"


#-------------------------------------------------
# READPOLYGON()
# Read Bezier control points from a file.
#
def ReadPolygon( filename ) :
    datafile = open(filename,'r');
    l = datafile.readline()
    degree = np.fromstring(l,sep=' ',dtype=int)
    BezierPts = np.fromfile(datafile,count=2*(degree+1),sep=' ',dtype=float)
    BezierPts = BezierPts.reshape(-1,2)
    return BezierPts


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
    #########
    ## TODO : Implement the De Casteljau algorithm.
    points = np.copy(BezierPts)
    for depth in range(0,k):
        for column in range(i,k-depth):
            points[column] = (1-t)*points[column] + t*points[column+1]
        # plt.plot( points[:,0], points[:,1], 'bo-' )
    return points[i]
    #########

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
    
    # degree of the curve (one less than the number of control points)
    degree = BezierPts.shape[0]-1
    
    # initialize curvepoints as zeros
    CurvePts = np.zeros([N,2])
    
    #########
    ## TODO : Compute N curve points for t varying uniformly in [0.0,1.0]
    #########
    #
    # hint1:
    # to generate the uniform sampling of the interval [0.0,1.0] with N elements, use:
    # >> samples = np.linspace(0.0,1.0,num=N)
    #
    samples = np.linspace(0.0,1.0,num=N)
    #
    # hint2:
    # to access and change i-th row of the matrix CurvePts, use:
    # >> CurvePts[i,:]
    #
    for i in range(0,N):
    	CurvePts[i,:] = DeCasteljau(BezierPts,degree,0,samples[i])
    # plt.plot( BezierPts[:,0], BezierPts[:,1], 'bo-' )
    #
    # hint3:
    # the actual point b_0^degree on the curve is computed by calling DeCasteljau
    # with arguments k=degree, i=0.
    #
    
    return CurvePts

def Polygons(BezierPts) :
    N = BezierPts.shape[0]-1
    # samples = np.linspace(0.0,1.0,num=N)
    t = 0.5
    points = np.copy(BezierPts)

    for depth in range(0,N):
        plt.plot(points[:,0], points[:,1], 'bo-', linewidth=(depth % 2.0)+0.1)

        for column in range(0,N-depth):
            points[column] = (1-t)*points[column] + t*points[column+1]
        points = points[:-1]
        # plt.plot( points[:,0], points[:,1], 'bo-' )
    plt.plot(points[:,0], points[:,1], 'bo-', linewidth=(depth % 2.0)+0.2)
#-------------------------------------------------
if __name__ == "__main__":
    
    # arg 1 : data name 
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
    else :
        dataname = "spiral" # simple, infinity, spiral, tuple

    # arg 2 : sampling density
    if len(sys.argv) > 2 :
        density = int(sys.argv[2])
    else :
        density = 50

    # filename
    filename = DATADIR + dataname + ".bcv"
    
    # check if valid datafile
    if not os.path.isfile(filename) :
        print ("error:  invalid dataname '" + dataname + "'")
        print ("usage:  python tp1.py  [simple,infinity,spiral,tuple]  [sampling_density]")
        
    else :    


        # plot
        fig = plt.gcf()
        fig.canvas.set_window_title('TP1 Bezier curves')
        plt.title(dataname+', '+str(density)+" pts")
        
        # set axes with equal proportions
        plt.axis('equal')
        
        # read control points
        BezierPts = ReadPolygon(filename)

        # compute curve points
        CurvePts = BezierCurve(BezierPts,density)

        # plot the control polygon
        # plt.plot( BezierPts[:,0], BezierPts[:,1], 'bo-' )
        
        # plot the curve
        plt.plot( CurvePts[:,0], CurvePts[:,1], 'r-' )
        
        #########
        ## TODO : Uncomment if you want to save the render as png image in the data/ dir
        #########
        #plt.savefig( DATADIR + dataname + ".png" )
        
        #########
        ## TODO : Compute intermediate polygons b_i^k for k=1,...,degree-1 and i=0,...,degree-k
        Polygons(BezierPts)
        # samples = np.linspace(0.0,1.0,num=density)
        # for i in range(0,density):
        #     points = np.copy(BezierPts)
        #     points[i,:] = DeCasteljau(points,points.shape[0]-1,0,samples[i])
        # plt.plot( BezierPts[:,0], BezierPts[:,1], 'bo-' )
        # for inter in range(1,density-1):
        # 	points = np.copy(BezierPts)
	       # 	DeCasteljau(points)
	       #  plt.plot( points[:,0], points[:,1], 'bo-' )
        #########
        
        #########
        ## TODO : Add plt.plot commands to plot the intermediate polygons
        #########
        
        plt.show()
