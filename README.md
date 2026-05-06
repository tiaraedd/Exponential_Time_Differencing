# Exponential Time Differencing (ETD) for Stiff PDEs

This project explores the mathematical background and implements a numerical framework for solving Partial Differential Equations (PDEs) that suffer from "stiffness". By utilizing **Exponential Time Differencing (ETD)**, the solver bypasses the stability constraints that often force traditional explicit schemes to use impractically small time steps.

## Key Features
* **Linear-Nonlinear Splitting**: The algorithm separates PDEs into stiff linear and non-stiff nonlinear components.
* **Exact Linear Integration**: Solves the linear portion of the system exactly via an integrating factor approach.
* **ETD-RK4 Integration**: Employs a fourth-order Runge-Kutta method to approximate the nonlinear integral term.
* **Spectral Acceleration**: Includes support for **Chebyshev nodes**, achieving exponential convergence and higher accuracy with fewer grid points compared to evenly spaced grids.
* **Generalised Boundary Handling**: Implements mechanisms for homogeneous and non-homogeneous Dirichlet conditions, as well as Neumann boundaries using ghost-point substitutions.

## Results
The provided implementation of the ETD method was tested on various PDE's, including the heat equation, Allen-Cahn equation, and Kuramoto-Sivashinsky equation. The results, as well as analyses of numerical convergence, computation speed, and comparison with other numerical solvers (including Crank-Nicolson and ADI), are summarized in the project write-up.

<!-- ## Mathematical Framework
The solver is designed for systems governed by the equation:
[cite_start]$$u_t = Lu + N(u, t)$$ [cite: 62]

[cite_start]Where $L$ is a linear operator and $N$ is the nonlinear portion[cite: 60]. The discrete update rule is formulated as:
[cite_start]$$u^{(m+1)} = e^{L\Delta t}u^{(m)} + e^{L\Delta t}\int_{0}^{\Delta t}e^{-Lt}N(u(t_m+\tau),t_m+\tau)d\tau$$ [cite: 64]

[cite_start]Matrix exponentials are computed using `scipy.linalg.expm`, which we found to be the most efficient approach for standard operators[cite: 75, 76]. -->

<!-- ## Implemented Models
[cite_start]The repository includes verified numerical solutions for the following systems[cite: 8, 201]:
1. [cite_start]**2D Heat Equation**: Used as a baseline for performance against Crank-Nicolson and ADI methods[cite: 82, 83, 84].
2. [cite_start]**Allen-Cahn Equation**: Benchmarked for phase separation modeling[cite: 169, 201].
3. [cite_start]**Fisher-KPP Equation**: For reaction-diffusion dynamics[cite: 201].
4. [cite_start]**Viscous Burgers' Equation**: Demonstrating stability in the presence of nonlinear convection[cite: 201].
5. [cite_start]**Kuramoto-Sivashinsky Equation**: Solved with periodic boundaries using Fast Fourier Transforms (FFT)[cite: 257, 258]. -->

<!-- ## Performance & Results
* [cite_start]**Stability**: Unlike standard explicit methods that "blow up" quickly, our ETD implementation remains stable and accurate for much larger time steps[cite: 96, 97].
* [cite_start]**Efficiency**: ETD outperforms traditional implicit solvers for coarser meshes where the initial cost of the matrix exponential solve is offset by the speed of the iterative steps[cite: 137, 139, 141].
* [cite_start]**Convergence**: Using Chebyshev grid points allowed us to achieve the same order of accuracy significantly faster than evenly spaced methods[cite: 286]. -->

## Authors
* **McKayla Davis** 
* **Tiara Eddington**

## References
*  Awad H. Al-Mohy and Nicholas J. Higham. Computing the action of the matrix exponential with
an application to exponential integrators. SIAM Journal on Scientific Computing, 33(2):488–511,
2011.
*  S. M. Cox and P. C. Matthews. Exponential time differencing for stiff systems. Journal of Computational Physics, 176(2):430–455, 2002.
*  Marlis Hochbruck and Alexander Ostermann. Exponential integrators. Acta Numerica, 19:209–286,
2010.
*  Aly-Khan Kassam and Lloyd N. Trefethen. Fourth-order time-stepping for stiff pdes. SIAM Journal
on Scientific Computing, 26(5):1214–1233, 2005
* Lloyd N. Trefethen. Spectral Methods in MATLAB. SIAM, Philadelphia, PA, 2000.