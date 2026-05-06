# [cite_start]Exponential Time Differencing (ETD) for Stiff PDEs [cite: 1]

[cite_start]This repository implements a numerical framework for solving Partial Differential Equations (PDEs) that suffer from "stiffness"[cite: 3, 12]. [cite_start]By utilizing **Exponential Time Differencing (ETD)**, the solver bypasses the stability constraints that often force traditional explicit schemes to use impractically small time steps[cite: 6, 12].

## Key Features
* [cite_start]**Linear-Nonlinear Splitting**: The algorithm separates PDEs into stiff linear and non-stiff nonlinear components[cite: 7, 15, 60].
* [cite_start]**Exact Linear Integration**: Solves the linear portion of the system exactly via an integrating factor approach[cite: 15, 55].
* [cite_start]**ETD-RK4 Integration**: Employs a fourth-order Runge-Kutta method to approximate the nonlinear integral term[cite: 67, 169].
* [cite_start]**Spectral Acceleration**: Includes support for **Chebyshev nodes**, achieving exponential convergence and higher accuracy with fewer grid points compared to evenly spaced grids[cite: 164, 165, 186].
* [cite_start]**Generalised Boundary Handling**: Implements mechanisms for homogeneous and non-homogeneous Dirichlet conditions, as well as Neumann boundaries using ghost-point substitutions[cite: 143, 144, 196, 198].

## Mathematical Framework
The solver is designed for systems governed by the equation:
[cite_start]$$u_t = Lu + N(u, t)$$ [cite: 62]

[cite_start]Where $L$ is a linear operator and $N$ is the nonlinear portion[cite: 60]. The discrete update rule is formulated as:
[cite_start]$$u^{(m+1)} = e^{L\Delta t}u^{(m)} + e^{L\Delta t}\int_{0}^{\Delta t}e^{-Lt}N(u(t_m+\tau),t_m+\tau)d\tau$$ [cite: 64]

[cite_start]Matrix exponentials are computed using `scipy.linalg.expm`, which we found to be the most efficient approach for standard operators[cite: 75, 76].

## Implemented Models
[cite_start]The repository includes verified numerical solutions for the following systems[cite: 8, 201]:
1. [cite_start]**2D Heat Equation**: Used as a baseline for performance against Crank-Nicolson and ADI methods[cite: 82, 83, 84].
2. [cite_start]**Allen-Cahn Equation**: Benchmarked for phase separation modeling[cite: 169, 201].
3. [cite_start]**Fisher-KPP Equation**: For reaction-diffusion dynamics[cite: 201].
4. [cite_start]**Viscous Burgers' Equation**: Demonstrating stability in the presence of nonlinear convection[cite: 201].
5. [cite_start]**Kuramoto-Sivashinsky Equation**: Solved with periodic boundaries using Fast Fourier Transforms (FFT)[cite: 257, 258].

## Performance & Results
* [cite_start]**Stability**: Unlike standard explicit methods that "blow up" quickly, our ETD implementation remains stable and accurate for much larger time steps[cite: 96, 97].
* [cite_start]**Efficiency**: ETD outperforms traditional implicit solvers for coarser meshes where the initial cost of the matrix exponential solve is offset by the speed of the iterative steps[cite: 137, 139, 141].
* [cite_start]**Convergence**: Using Chebyshev grid points allowed us to achieve the same order of accuracy significantly faster than evenly spaced methods[cite: 286].

## Authors
* [cite_start]**McKayla Davis** [cite: 2, 20, 205]
* [cite_start]**Tiara Eddington** [cite: 2, 20, 205]

## References
* [cite_start][AMH11] Al-Mohy & Higham, "Computing the action of the matrix exponential"[cite: 297].
* [cite_start][CM02] Cox & Matthews, "Exponential time differencing for stiff systems"[cite: 299].
* [cite_start][HO10] Hochbruch & Ostermann, "Exponential integrators"[cite: 301].
* [cite_start][KT05] Kassam & Trefethen, "Fourth-order time-stepping for stiff pdes"[cite: 303].
* [cite_start][Tre00] Trefethen, "Spectral Methods in MATLAB"[cite: 305].