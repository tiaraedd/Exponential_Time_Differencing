import numpy as np 
from numpy import linalg as la

def cheb(N):
    ''' Given a specified number of points N, computes the Chebyshev differentiation matrix
        a la Trefethen's algorithm from his textbook "Spectral Methods in MATLAB".

        Parameters
            N (int) : The number of Chebyshev nodes to use in the discretization
        
        Returns
            D ((N+1), (N+1) ndarray) : the Chebyshev differentiation matrix 
            x ((N+1, ) ndarray) : the spatial mesh nodes (Chebyshev nodes)
    '''

    # explicit case for the trivial (1-point) mesh
    if N == 0:
        D = np.array(0)
        x = np.array(1)
        return D, x 

    # define points in mesh (chebyshev points from -1 to 1)
    lin = np.arange(0, N+1, dtype=np.float64)
    x = np.cos(np.pi * lin / N)       

    # compute a vector of coefficients for the derivative values
    # (weights with alternating signs)
    c = np.ones(N+1, dtype=np.float64)
    c[0], c[-1] = 2., 2.
    c *= (-1.)**lin

    # compute the distance matrix
    dX = x[:, None] - x[None, :]

    # compute differentiation matrix
    D = (c[:, None] / c) / (dX + np.eye(N+1))       # this is a neat way to compute an outer product
    D -= np.diag(np.sum(D,axis=1))
    
    # Here we return the differentiation matrix and the Chebyshev points,
    # numbered from x_0 = -1 to x_N = 1
    return D[::-1,::-1], x[::-1]


def etd_solve(L, u0, T, nt):
    ''' Use ETD to solve a PDE of the form u_t = L(u), where L is a linear operator (i.e. a matrix apprximation of
        the appropriate derivative operator). 
        ** IMPORTANT : This currently only works for homogenous Dirichley BC's. The function slices out the interior 
                       of the linear operator so as not to affect the boundaries
        
        Parameters:
            L (ndarray) : the linear operator defining/approximating the right side of the PDE
            u0 (ndarray) : the initial state of the system 
            T (float) : the final time value 
            nt (int) : the desired number of time steps
        
        Returns:
            U (ndarray) : the approximate solution -- each row corresponds to a time step 
            t_vals (ndarray) : the linspace of t-values 
    '''
    
    t_vals, dt = np.linspace(0, T, nt, retstep=True)    # generate t values 
    U = np.empty((nt, len(u0)))                         # initialize solution array 
    U[0] = u0                                           # set initial condition
    U[:, 0] = 0.                                        # fix left boundary
    U[:, -1] = 0.                                       # fix right boundary

    # approximate the matrix exponential e^(dt*L)
    # testing shows that this is fastest by diagonalizing the matrix
    L_int = L[1:-1, 1:-1]                               # slice out interior of operator
    lam, V = la.eig(dt * L_int)
    exp_mat = V @ np.diag(np.exp(lam)) @ la.inv(V)

    # solver loop 
    for i in range(1, nt):
        U[i, 1:-1] = exp_mat @ U[i-1, 1:-1]
    
    return U, t_vals
