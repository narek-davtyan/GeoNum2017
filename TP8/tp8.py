#-------------------------------------------------
#
#  TP8 : Uniform B-splines as Subdivision Surfaces
#  http://tiborstanko.sk/teaching/geo-num-2017/tp8.html
#  [31-Mar-2017]
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
# READPOINTS()
# Read control network (grid) of points from a file.
#
# Input
#     datafile  :  .net file to be read
#
# Output 
#     M         :  matrix, 3 x m x n, control grid (all three coordinates)
#     u_closed  :  True  = surface is closed in u-direction
#                  False = surface is open   in u-direction
#     v_closed  :  True  = surface is closed in v-direction
#                  False = surface is open   in v-direction
#
def ReadPoints( datafile ) :
    m,n,u_closed,v_closed = np.fromfile(datafile,count=4,sep=' ',dtype=int)
    M = np.fromfile(datafile,count=3*m*n,sep=' ',dtype=float)
    M = M.reshape(-1,3).transpose().reshape(3,m,n)
    return M, u_closed==1, v_closed==1

#-------------------------------------------------
# GENERATERANDOMTERRAIN()
# Generate random terrain and save it to data/terrain.net
#
# Input
#      N    :  grid size in X, Y
#    Amp    :  amplitude in Z
#
def GenerateRandomTerrain(N,Amp) :
    X,Y = np.meshgrid(range(N-1,-1,-1),range(N))
    V = np.array([X,Y,Amp*np.random.rand(N,N)]).transpose().reshape(-1,3)
    f = open(TP+"data/terrain.net","wa")
    f.write('%d %d %d %d \n' % (N,N,0,0))
    np.savetxt(f,V,fmt='%.4f')
    f.close()

#-------------------------------------------------
# SUBDIVIDE()
# Perform one subdivision step.
#
# Input
#        M0   :  3 x M x N, control network of points
#  u_closed   :  True if closed in u-direction, False otherwise
#  v_closed   :  True if closed in v-direction, False otherwise
#
# Output
#       M1    :  array with subdivided network
#
def Subdivide( M0, u_closed, v_closed ) :
    
    # dim is equal to 3
    dim, m, n = M0.shape
    
    # upsample
    M1 = np.zeros([dim, \
        2*m if u_closed else 2*(m-1), \
        2*n if v_closed else 2*(n-1)])
    
    # calculating and applying subdivision masks
    for i in range(0,m if u_closed else (m-1)) :
        for j in range(0,n if v_closed else (n-1)) :
            
            i_plus_one = (i+1)%m if u_closed else (i+1)
            j_plus_one = (j+1)%n if v_closed else (j+1)
            
            A = M0[:,i,j]
            B = M0[:,i_plus_one,j]
            C = M0[:,i,j_plus_one]
            D = M0[:,i_plus_one,j_plus_one]

            M1[:,2*i+0,2*j+0] = (1.0/16.0)*(9.0*A +3.0*B +3.0*C +1.0*D)
            M1[:,2*i+1,2*j+0] = (1.0/16.0)*(3.0*A +9.0*B +1.0*C +3.0*D)
            M1[:,2*i+0,2*j+1] = (1.0/16.0)*(3.0*A +1.0*B +9.0*C +3.0*D)
            M1[:,2*i+1,2*j+1] = (1.0/16.0)*(1.0*A +3.0*B +3.0*C +9.0*D)    

    return M1

def Sub_C_C( M0, u_closed, v_closed ) :

    # dim is equal to 3
    dim, m, n = M0.shape

    # generate face points
    M1 = np.zeros([dim, \
        m if u_closed else (m-1), \
        n if v_closed else (n-1)])
    for i in range(0,m if u_closed else (m-1)) :
        for j in range(0,n if v_closed else (n-1)) :

            i_plus_one = (i+1)%m if u_closed else (i+1)
            j_plus_one = (j+1)%n if v_closed else (j+1)
            
            A = M0[:,i,j]
            B = M0[:,i_plus_one,j]
            C = M0[:,i,j_plus_one]
            D = M0[:,i_plus_one,j_plus_one]

            M1[:,i+0,j+0] = (1.0/4.0)*(A + B + C + D)
    print(M0)
    print(M1)

    # generate edge points
    M2 = np.zeros([dim, \
        m/2+m%2 if u_closed else (m/2+m%2-1), \
        n/2+n%2 if v_closed else (n/2+n%2-1)])
    M3 = np.zeros([dim, \
        n/2+m%2 if v_closed else (n/2+n%2-1), \
        m/2+n%2 if u_closed else (m/2+m%2-1)])
    for i in range(0,m-1 if u_closed else (m-2)) :
        for j in range(0,n-1 if v_closed else (n-2)) :
            
            i_plus_one = (i+1)%m if u_closed else (i+1)
            j_plus_one = (j+1)%n if v_closed else (j+1)
            i_plus_two = (i+2)%m if u_closed else (i+2)
            j_plus_two = (j+2)%n if v_closed else (j+2)
            
            A = M0[:,i,j]
            B = M0[:,i_plus_one,j]
            C = M0[:,i,j_plus_one]
            D = M0[:,i_plus_one,j_plus_one]
            E = M0[:,i,j_plus_two]
            F = M0[:,i_plus_one,j_plus_two]
            G = M0[:,i_plus_two,j]
            H = M0[:,i_plus_two,j_plus_one]

            M2[:,i/2+i%2,j/2+j%2] = (1.0/4.0)*(A + B + 6.0*C + 6.0*D + E + F)
            M3[:,i/2+i%2,j/2+j%2] = (1.0/4.0)*(A + C + 6.0*B + 6.0*D + G + H)

    print(M2)
    print(M3)

    # modify the existing vertices
    # (Q/n) + (2R/n) + (S(n-3)/n) , where valence n=4, hence (Q/4) + (R/2) + (S/4),
    # with Q - the average of the surrounding face points, 
    # R - the average of all surround edge midpoints, and 
    # S - the control point
    M4 = np.zeros([dim, \
        n if v_closed else (n-1), \
        m if u_closed else (m-1)])
    for i in range(0,m-1 if u_closed else (m-2)) :
        for j in range(0,n-1 if v_closed else (n-2)) :
            
            i_plus_one = (i+1)%m if u_closed else (i+1)
            j_plus_one = (j+1)%n if v_closed else (j+1)

            A = M1[:,i,j]
            B = M1[:,i_plus_one,j]
            C = M1[:,i,j_plus_one]
            D = M1[:,i_plus_one,j_plus_one]
            Q = (1.0/4.0)*(A + B + C + D)

            i_plus_one = (i+1)%(m/2+m%2) if u_closed else (i+1)
            j_plus_one = (j+1)%(n/2+n%2) if v_closed else (j+1)

            E = M2[:,i,j]
            F = M2[:,i_plus_one,j]
            G = M3[:,i,j]
            H = M3[:,i,j_plus_one]
            R = (1.0/4.0)*(E + F + G + H)
            S = M0[:,i_plus_one,j_plus_one]

            M4[:,i,j] = (1.0/4.0)*Q + (1.0/2.0)*R + (1.0/4.0)*S
            
#-------------------------------------------------
if __name__ == "__main__":
    
    # arg 1 : data name
    if len(sys.argv) > 1 :
        dataname = sys.argv[1]
        if dataname == "terrain" :
            # First parameter modifies the length/width
            # Second parameter modifies the height
            GenerateRandomTerrain(59,1.0)
            #GenerateRandomTerrain(20,3.0)
    else :
        # torus,cylinder,grid,terrain
        dataname = "torus"

    # arg 2 : subdivision depth
    if len(sys.argv) > 2 :
        # subdivision_depth=1
        depth = int(sys.argv[2])
    else :
        depth = 1

    # arg 3 : algo choice
    if len(sys.argv) > 3 :
        # algo_cc=False
        algo_cc = True
    else :
        algo_cc = False

    # filename
    filename = TP+"data/"+dataname+".net"
    
    # check if valid datafile
    if not os.path.isfile(filename) :
        print " error   :  invalid dataname '" + dataname + "'"
        print " usage   :  tp8.py  [torus,cylinder,grid,terrain]  [subdivision_depth=1]"
        print " example :  python tp8.py torus 3"
        
    else :
        # init Viewer
        viewer = Viewer("TP8 : Uniform B-spline SubSurf["+dataname+"]",[1200,800])
        
        # open the datafile
        datafile = open(filename,'r');
        
        # read control points and u_closed, v_closed
        M, u_closed, v_closed = ReadPoints( datafile )
        
        # add wireframe : control net
        viewer.add_patch( M[0,:,:], M[1,:,:], M[2,:,:], wireframe=True)
                
        # iterative subdivision
        for d in range(depth) :
            Sub_C_C(M,u_closed,v_closed)
            M = Subdivide(M,u_closed,v_closed) if not algo_cc else Sub_C_C(M,u_closed,v_closed)
        
        # u_closed : make rows periodic
        if u_closed :
            rows = np.append( np.arange(M.shape[1]), 0)
            M = M[:,rows,:]
            
        # v_closed : make cols periodic
        if v_closed :
            cols = np.append( np.arange(M.shape[2]), 0)
            M = M[:,:,cols]
        
        # add patch : subdivision surface
        viewer.add_patch( M[0,:,:], M[1,:,:], M[2,:,:])
        
        # display the viewer
        viewer.render()
