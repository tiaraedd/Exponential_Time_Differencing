import numpy as np
from matplotlib import pyplot as plt
from matplotlib.animation import FuncAnimation


def crank_heat_equation(a, h_x, h_t, U_0, x_0, x_n, t_0, t_n, filename = 'heat_solution.mp4'):
    '''
    Solves the 1-d heat equation with dirichlet boundary conditions using the crank-nicolson scheme.
    INPUTS:
    a (float): a from the heat equation
    h_x (float): the step size for space
    h_t (float): the step size for time
    U_0 (func): Function representing initial condition
    
    RETURNS:
    U (np.array): Solution for u at each timestep
    mu (float): CFL value
    
    '''
    # Create the grid
    x = np.linspace(x_0, x_n, int((x_n - x_0)/h_x) + 1)
    T = np.linspace(t_0, t_n, int((t_n - t_0)/h_t) + 1)
    
    # Build P and Q matrices for fDM method
    mu = a * h_t/(h_x**2)
    P = np.diag(mu/2 * np.ones(len(x) - 1), k = -1) + np.diag(np.ones_like(x) - mu, k = 0) + np.diag(mu/2 * np.ones(len(x) - 1), k = 1)
    Q = np.diag(-mu/2 * np.ones(len(x) - 1), k = -1) + np.diag(np.ones_like(x) + mu, k = 0) + np.diag(-mu/2 * np.ones(len(x) - 1), k = 1)
    
    # Fix boundary Conditions
    P[0], P[-1] = 0, 0
    P[0,0], P[-1,-1] = 1, 1
    Q[0], Q[-1] = 0, 0
    Q[0,0], Q[-1,-1] = 1, 1
    
    # Iterate in time
    U = [U_0(x)]
    for n in range(len(T)-1):
        left = P @ U[-1]
        
        # Fix boundary conditions in time
        left[0], left[-1] = 0, 0
            
        solution = np.linalg.solve(Q, left)
        U.append(solution)
    U = np.array(U)
        
    # Print the CFL relationship
    print(f"CFL condition: a * h_t/ (h_x**2) = {mu}")
    
    # Create animation of solution
    # Create a figure and axis object
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.xlim([x_0, x_n])
    plt.ylim([np.min(U) - 1, np.max(U) + 1])
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')

    # Create an empty line object
    line, = ax.plot([], [], 'r-')
    plt.suptitle(f"Crank-Nicolson Solution for heat equation with CFL {mu}")

    # Make update function
    def update(t):
        line.set_data(x, U[t])
        ax.set_title(f"t = {T[t]}")

    # Solve
    ani = FuncAnimation(fig, update, frames=np.arange(0,len(T)), interval=(4000*(T[-1]-T[0]))/(len(T)-1))
    ani.save(filename)
    plt.close()
        
    # Return solution and CFL condition
    return U, mu


def crank_wave_equation(a, h_x, h_t, U_0, x_0, x_n, t_0, t_n, filename = 'wave_solution.mp4'):
    '''
    Solves the 1 way wave equation with dirichlet boundary conditions using the crank-nicolson scheme.
    INPUTS:
    a (float): a from the heat equation
    h_x (float): the step size for space
    h_t (float): the step size for time
    U_0 (func): Function representing initial condition
    
    RETURNS:
    U (np.array): Solution for u at each timestep
    mu (float): CFL value
    
    '''
    # Create the grid
    x = np.linspace(x_0, x_n, int((x_n - x_0)/h_x) + 1)
    T = np.linspace(t_0, t_n, int((t_n - t_0)/h_t) + 1)
    
    # Build P and Q matrices for fDM method
    mu = a * h_t/(h_x)
    P = np.diag(mu/4 * np.ones(len(x) - 1), k = -1) + np.diag(np.ones_like(x), k = 0) + np.diag(-mu/4 * np.ones(len(x) - 1), k = 1)
    Q = np.diag(-mu/4 * np.ones(len(x) - 1), k = -1) + np.diag(np.ones_like(x), k = 0) + np.diag(mu/4 * np.ones(len(x) - 1), k = 1)
    
    # Fix boundary Conditions
    P[0], P[-1] = 0, 0
    P[0,0], P[-1,-1] = 1, 1
    Q[0], Q[-1] = 0, 0
    Q[0,0], Q[-1,-1] = 1, 1
    
    # Iterate in time
    U = [U_0(x)]
    for n in range(len(T)-1):
        left = P @ U[-1]
        
        # Fix boundary conditions in time
        left[0], left[-1] = 0, 0
            
        solution = np.linalg.solve(Q, left)
        U.append(solution)
    U = np.array(U)
        
    # Print the CFL relationship
    print(f"CFL condition: a * h_t/ h_x = {mu}")
    
    # Create animation of solution
    # Create a figure and axis object
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    plt.xlim([x_0, x_n])
    plt.ylim([np.min(U) - 1, np.max(U) + 1])
    ax.set_xlabel('x')
    ax.set_ylabel('u(x,t)')

    # Create an empty line object
    line, = ax.plot([], [], 'r-')
    plt.suptitle(f"Crank-Nicolson Solution for one way wave equation with CFL {mu}")

    # Make update function
    def update(t):
        line.set_data(x, U[t])
        ax.set_title(f"t = {T[t]}")

    # Solve
    ani = FuncAnimation(fig, update, frames=np.arange(0,len(T)), interval=(4000*(T[-1]-T[0]))/(len(T)-1))
    ani.save(filename)
    plt.close()
        
    # Return solution and CFL condition
    return U, mu