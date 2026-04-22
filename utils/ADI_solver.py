import numpy as np
from matplotlib import pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML
from matplotlib.animation import FuncAnimation
import time

def ADI_Solver(alpha, h_x, h_y, h_t, plot = False, filename = 'adi_solution.mp4', initial = None, truth = None):
    '''
    Solves the 2-d heat equation with homogeneous boundary conditions using the alternating direction implicit method.
    
    INPUTS:
    alpha (float): a from the heat equation
    h_x (float): the step size for x direction
    h_x (float): the step size for y direction
    h_t (float): the step size for time
    plot(bool): whether or not you would like an animation drawn for your solution
    filename (str): the filename for the animation
    initial (func): a function taking in X, Y representing the initial condition of the solution
    truth (func): a function taking in X,Y representing the true solution at the final time step
    
    RETURNS:
    U (np.array): solution for u at each timestep
    norm (float): the error
    time (float): the time it took to solve the system
    
    '''
    # Initialize intial condition if needed
    if initial == None:
        initial = lambda x,y: 1 * np.sin(np.pi * x) * 1 * np.sin(np.pi * y)
        truth = lambda X, Y, t: np.exp(-2 * np.pi**2 * alpha * t) * np.sin(np.pi * X) * np.sin(np.pi * Y)
    
    # Define mesh
    x = np.linspace(0,1,int(1/h_x) + 1)
    y = np.linspace(0,1,int(1/h_y) + 1)
    t_vals = np.linspace(0,1,int(1/h_t) + 1)

    X,Y = np.meshgrid(x,y)

    # Derive ADI matrices
    A_1 = alpha/(h_x**2) * (np.diag(np.ones_like(x)[1:], k = -1) + np.diag(np.ones_like(x) * -2, k = 0) + np.diag(np.ones_like(x)[1:], k = 1))
    A_2 = alpha/(h_y**2) * (np.diag(np.ones_like(y)[1:], k = -1) + np.diag(np.ones_like(y) * -2, k = 0) + np.diag(np.ones_like(y)[1:], k = 1))

    # Define boundary conditions (not needed for homogeneous case)
    # B_top1 = (np.diag(1 * np.ones_like(y)[1:]/h_x, k = -1) + np.diag(np.ones_like(y) * 2 * (1 - 1/h_x), k = 0) + np.diag(np.ones_like(y)[1:]/h_x, k = 1))[1:-1,1:-1]/4
    # B_top2 = (np.diag(-1 * np.ones_like(y)[1:]/h_x, k = -1) + np.diag(np.ones_like(y) * 2 * (1 + 1/h_x), k = 0) + np.diag(-np.ones_like(y)[1:]/h_x, k = 1))[1:-1,1:-1]/4

    U_pred = [initial(X,Y)]
    start = time.time()
    
    # Iterate through each time step and solve the system
    for t in t_vals[1:]:
        # Evaluate the half step
        v_half = np.ones_like(U_pred[-1])
        v_half[1:-1,1:-1] = np.linalg.solve(np.eye(len(A_1)-2) - h_t/2 * A_1[1:-1, 1:-1], ((np.eye(len(A_2) - 2) + h_t/2 * A_2[1:-1,1:-1]) @ U_pred[-1][1:-1,1:-1]).T).T

        # Adjust boundary condition
        # v_half[0,1:-1] = B_top1 @ u_true(t - h_t, X[0,1:-1], Y[0,1:-1]) + B_top2 @ u_true(t, X[0,1:-1], Y[0,1:-1])
        # v_half[-1,1:-1] = B_top1 @ u_true(t - h_t, X[-1,1:-1], Y[-1,1:-1]) + B_top2 @ u_true(t, X[-1,1:-1], Y[-1,1:-1])
        v_half[0,1:-1] = 0
        v_half[-1,1:-1] = 0

        # Solve the system at the time step
        v = np.ones_like(U_pred[-1])
        v[1:-1,1:-1] = np.linalg.solve(np.eye(len(A_2) - 2) - h_t/2 * A_2[1:-1, 1:-1], ((np.eye(len(A_1) - 2) + h_t/2 * A_1[1:-1, 1:-1]) @ v_half[1:-1,1:-1].T).T)

        # Set boundary conditions
        v[0,:] = 0
        v[-1,:] = 0
        v[1:-1,0] = 0
        v[1:-1,-1] = 0

        # Save time step
        U_pred.append(v)
    
    end_time = time.time() - start

    
    if plot == True:
        # Create 3D plot
        fig = plt.figure(figsize=(8,6))
        ax1 = fig.add_subplot(111, projection='3d')

        ax1.set_zlim3d([-1.1,1.1])

        # Create empty line objects
        pred = ax1.plot_surface(X, Y, U_pred[0], cmap = 'inferno')

        ax1.set_title('Predicted Solution')
        # Make update function to updata all line plots
        def update(t):
            global true, pred

            pred.remove()

            pred = ax1.plot_surface(X, Y, U_pred[t], cmap='inferno', vmax = np.max(U_pred), vmin = np.min(U_pred))

            return pred


        # Make figure
        ani = FuncAnimation(fig, update, frames=range(len(t_vals)), interval=40)
        ani.save(filename)

        plt.close()

    return U_pred, np.linalg.norm(U_pred[-1] - truth(X,Y,t)), end_time


