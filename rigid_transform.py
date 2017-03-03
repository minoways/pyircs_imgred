import numpy as np
import math 

def rigid_transform(A, B):

    assert len(A) == len(B)

    N = A.shape[0]; # total points

    if N > 2:

        centroid_A = np.mean(A, axis=0)
        centroid_B = np.mean(B, axis=0)
    
        # centre the points
        AA = A - np.tile(centroid_A, (N, 1))
        BB = B - np.tile(centroid_B, (N, 1))

        # dot is matrix multiplication for array
        H = np.transpose(AA) * BB

        U, S, Vt = np.linalg.svd(H)

        R = Vt.T * U.T

        t = -R*centroid_A.T + centroid_B.T

    elif N == 2:
        
        dA = A[1] - A[0]
        dB = B[1] - B[0]

        a = (dA[0,0]*dB[0,0] + dA[0,1]*dB[0,1]) / (dA[0,0]**2 + dA[0,1]**2)
        b = (dA[0,0]*dB[0,1] - dB[0,0]*dA[0,1]) / (dA[0,0]**2 + dA[0,1]**2)
        m = '%f %f; %f %f' % (a, -b, b, a)
        
        R = np.matrix(m)
        t = np.mean((B.T - R*A.T), axis=1)
    
    else:

        R = np.matrix('1.0 0; 0 1.0')
        t = (B[0] - A[0]).T

    return R, t
