##################################################
######         DAVTYAN && HADJADJ           ######
##################################################
#-------------------------------------------------
#
#  TP7 : B-spline surfaces
#  http://tiborstanko.sk/teaching/geo-num-2017/tp7.html
#  [24-Mar-2017]
#
#-------------------------------------------------
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
#-------------------------------------------------

import sys, os
import numpy as np
from viewer import Viewer

TP = os.path.dirname(os.path.realpath(__file__)) + "/"

#-------------------------------------------------
# READBSPLINEMESHWITHKNOTS()
# Read B-spline surface data from a .bspline of a .nurbs file.
#
# Input
#     datafile  :  file to be read
#
# Output 
#     M         :  matrix, 3 (or 4) x (m+1) x (n+1), control net (all three (four) coordinates)
#     U         :  knot sequence in the u-direction
#     V         :  knot sequence in the v-direction
#
def ReadBSplineMeshWithKnots( datafile, nurbs=False ) :
    # if nurbs, read four coordinates; otherwise read three
    if nurbs : 
        dim=4
    else : 
        dim=3
    m,n,k,l = np.fromfile(datafile,count=4,sep=' ',dtype=int)
    M = np.fromfile(datafile,count=dim*m*n,sep=' ',dtype=float)
    M = M.reshape(-1,dim).transpose().reshape(dim,m,n)
    U = np.fromfile(datafile,count=k,sep=' ',dtype=float)
    V = np.fromfile(datafile,count=l,sep=' ',dtype=float)
    return M, U, V

#-------------------------------------------------
# DEBOOR( ... )
# Perform the De Boor's algorithm recursively.
#
# Input
#    ControlPts :  (n+1) x 2 matrix of control points
#    Knots      :  (m+1) x 1 vector of knots
#    r          :  upper index of the computed point (depth of the algorithm)
#    j          :  lower index of the computed point
#    t          :  curve parameter in [t_i,t_i+1]
#
# Output
#   Point d_j^r from the De Boor algorithm.
#

def DeBoor( ControlPts, Knots, r, j, t ) :
    k = Knots.shape[0] - ControlPts.shape[0] -1
    if r==0 :
        return ControlPts[j]
    else :
        # (1-w_j_k-r+1(t)) * d_j-1^(r-1) 
        #  + w_j_k-r+1(t) * d_j^(r-1) 
        return (1-ComputeW(Knots, j, k-r+1, t )) * \
        DeBoor( ControlPts, Knots, r-1, j-1, t )  +  \
        ComputeW(Knots, j, k-r+1, t ) * \
        DeBoor( ControlPts, Knots, r-1, j, t )

def ComputeW( Knots, i, k, t) :
    if Knots[i] < Knots[i+k] :
        return (t - Knots[i])/(Knots[i+k] - Knots[i])
    else:
        return 0

#-------------------------------------------------
# DEBOORSURF( ... )
# Recursive De Boor's algorithm for surfaces.
#
# Input
#    M     :  (m+1) x (n+1) coordinate matrix, control points net
#    U     :  vector of knots in the direction u
#    V     :  vector of knots in the direction v
#    r, s  :  upper indices of the computed point (depth of the algorithm)
#    i, j  :  lower indices of the computed point
#    u, v  :  parameters, [u,v] in [ U_i, U_i+1 ] x [ V_j, V_j+1 ]
#
# Output
#   Point d_(i,j)^(r,s) from De Boor's algorithm.
# 
def DeBoorSurf( M, U, V, r, s, i, j, u, v ) :
    m = M.shape[0]-1
    n = M.shape[1]-1

    tmp = np.zeros([n+1])
    k=0
    # Performing DeBoor n+1 times in direction u
    for d in range(0,n+1) :
        tmp[k] = DeBoor(M[:,k], U, r, i, u)
        k+=1
    # Performing DeBoor 1 time in direction v
    return DeBoor(tmp, V, s, j, v)

#-------------------------------------------------
if __name__ == "__main__":
    
    # arg 1 : data name
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
    else :
        # simple (bspline), torus (bspline+nurbs), hemi (nurbs)
        dataname = "simple" 

    # arg 2 : sampling density
    if len(sys.argv) > 2 :
        samples = int(sys.argv[2])
    else :
        samples = 11
        
    # arg 3 : nurbs
    if len(sys.argv) > 3 :
        nurbs = True
        ext = "nurbs"
    else :
        nurbs = False
        ext = "bspline"

    # filename
    filename = TP+"data/"+dataname+"."+ext
    
    # check if valid datafile
    if not os.path.isfile(filename) :
        print (" error   :  invalid dataname '" + dataname + "'")
        print (" usage   :  tp7.py  [simple,torus]  [sampling_density]")
        print (" example :  python tp7.py torus 20")
        print (" usage   :  tp7.py  [torus, hemi]  [sampling_density] nurbs")
        print (" example :  python tp7.py torus 22 nurbs")
        
        
    else :
        # init Viewer
        viewer = Viewer("TP7 : B-spline surfaces ["+dataname+"]",[1200,800])
        
        # open the datafile
        datafile = open(filename,'r');
        
        # read control points and knot sequences
        M, U, V = ReadBSplineMeshWithKnots( datafile, nurbs )

        # coordinate matrices
        Mx = M[0,:,:]
        My = M[1,:,:]
        Mz = M[2,:,:]
        
        # add wireframe : control net
        viewer.add_patch( Mx, My, Mz, wireframe=True)
        
        # NURBS weights : Multiply Mx, My, Mz by Mw
        if nurbs :
            Mw = M[3,:,:]
            Mx = np.multiply(Mx,Mw)
            My = np.multiply(My,Mw)
            Mz = np.multiply(Mz,Mw)
        
        # add control net wireframe to the viewer
        viewer.add_patch( Mx, My, Mz, wireframe=True)

        m = Mx.shape[0]-1    # m+1 points in u-direction
        n = Mx.shape[1]-1    # n+1 points in v-direction
        
        k = U.shape[0]-1     # k+1 knots in u-direction
        l = V.shape[0]-1     # l+1 knots in v-direction
        
        du = k-m-1           # degree in u-direction
        dv = l-n-1           # degree in v-direction
        
        # loop over segments in u-direction
        for i in range(du,k-du) :
            
            # check if the segment is non-degenerate
            if U[i] == U[i+1] :
                continue
            
            # loop over segments in v-direction
            for j in range(dv,l-dv) :
                
                # check if the segment is non-degenerate
                if V[j] == V[j+1] :
                    continue
            
                ##   we are now evaluating a surface patch defined over 
                ##   [ U_i, U_i+1 ] x [ V_j, V_j+1 ]
                
                # initialize patch points : three coordinate matrices
                Sx = np.zeros([samples,samples])
                Sy = np.zeros([samples,samples])
                Sz = np.zeros([samples,samples])
                if nurbs :
                    Sw = np.zeros([samples,samples])

                ## Compute patch points using DeBoorSurf.
                ##   Use a double loop with uniform sampling
                ##    - loop over u in np.linspace(U_i, U_i+1,num=samples)
                ##    -- loop over v in np.linspace(V_j, V_j+1,num=samples)

                si=0
                for u in np.linspace(U[i],U[i+1],num=samples) :
                    sj=0
                    for v in np.linspace(V[j],V[j+1],num=samples) :
                        Sx[si,sj] = DeBoorSurf(Mx,U,V,du,dv,i,j,u,v)
                        Sy[si,sj] = DeBoorSurf(My,U,V,du,dv,i,j,u,v)
                        Sz[si,sj] = DeBoorSurf(Mz,U,V,du,dv,i,j,u,v)
                        if nurbs :
                            Sw[si,sj] = DeBoorSurf(Mw,U,V,du,dv,i,j,u,v)
                        sj+=1   
                    si+=1

                
                ## NURBS Weights : Divide Sx, Sy, Sz by Sw
                if nurbs :
                    Sx = np.divide(Sx,Sw)
                    Sy = np.divide(Sy,Sw)
                    Sz = np.divide(Sz,Sw)
                
                # after the Sx, Sy, Sz have been calculated :
                # add current patch to the viewer
                viewer.add_patch(Sx,Sy,Sz)
                
            # END for j in range( degreeV, n-degreeV )
        # END for i in range( degreeU, m-degreeU )
        
        # display the viewer
        viewer.render()
