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
import sys

TP = os.path.dirname(os.path.realpath(__file__)) + "/"
DATADIR = filename = TP+"data/"

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
def DeCasteljauIterative( BezierPts, k, i, t ) :

    # compute deeper level points and write them over the old values
    # b_column = (1-t)*b_column + t*b_(column-1)
    for depth in range(0,k):
        for column in range(i,k-depth):
            BezierPts[column] = (1-t)*BezierPts[column] + t*BezierPts[column+1]
    
    # return the i-th point
    return BezierPts[i]
def DeCasteljauRecursive( BezierPts, k, i, t ) :

    if k!=0 :
        # return the i-th point of intermediary levels
        return (1-t)*DeCasteljauRecursive(BezierPts, k-1, i, t) \
         + t*DeCasteljauRecursive(BezierPts, k-1, i+1, t)
    else :
        # return the i-th point of first levels
        return BezierPts[i]

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

    # generate the uniform sampling of the interval [0.0,1.0] with N elements
    samples = np.linspace(0.0,1.0,num=N)
    # print (samples)

    # compute N curve points for t varying uniformly in [0.0,1.0]
    for i in range(0,N):
        points = np.copy(BezierPts)
        if rec == True :
            CurvePts[i] = DeCasteljauRecursive(points,degree,0,samples[i])
        else :
            CurvePts[i] = DeCasteljauIterative(points,degree,0,samples[i])
    
    return CurvePts

#-------------------------------------------------
# POLYGONPRINTER( ... )
# Compute intermediate polygons and
# plot the corresponding segments
# b_i^k for k=1,...,degree-1 and i=0,...,degree-k
#
# Input
#    BezierPts :  (degree+1) x 2 matrix of Bezier control points
#
def PolygonPrinter( BezierPts ) :
    # get number of segments
    N = BezierPts.shape[0]-1

    # fix the value of t
    t = 0.5

    # copy the original array with Bezier points
    points = np.copy(BezierPts)

    for depth in range(0,N):
        
        # compute points of the next polygon
        for column in range(0,N-depth):
            points[column] = (1-t)*points[column] + t*points[column+1]

        # omit the last point since it's no longer useful
        points = points[:-1]
        
        # print all intermediate polygons
        plt.plot(points[:,0], points[:,1], '-o', 
            linewidth=1.2 if depth%2.0==0 else 0.6)

#-------------------------------------------------
if __name__ == "__main__":
    
    # arg 1 : data name 
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
    else :
        dataname = "infinity" # simple, infinity, spiral, tuple

    # arg 2 : sampling density
    if len(sys.argv) > 2 :
        density = int(sys.argv[2])
    else :
        density = 50

    # arg 3 : recursive or iterative DeCasteljau
    if len(sys.argv) > 3 :
        rec = (True if sys.argv[3]=="rec" else False)
    else :
        rec = True

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

        # print the control polygon
        plt.plot(BezierPts[:,0], BezierPts[:,1], '-o', linewidth=1.5)

        # plot the curve
        plt.plot( CurvePts[:,0], CurvePts[:,1], 'r-' )
        
        # save the render as png image in the data/ dir
        #plt.savefig( DATADIR + dataname + ".png" )
        
        # compute and plot intermediate polygons 
        # b_i^k for k=1,...,degree-1 and i=0,...,degree-k
        PolygonPrinter(BezierPts)

        plt.show()
