# The book: *“Parallel Software for the Robust Multigrid Technique”* by S.&nbsp;I.&nbsp;Martynenko

The repository contains examples and programs from the book *“Parallel Software for the Robust Multigrid Technique”* by S.&nbsp;I.&nbsp;Martynenko, Triumph, Moscow, 2021. Complete text of the book in Russian is available via the link. Abstract in English, theoretical background and description of the source code files are provided below.

## Abstract of the book

This book covers parallel software for the Robust Multigrid Technique and presents basic concepts for black box software development in the computational continuum mechanics (thermal conductivity, chemical hydrodynamics, convective heat transfer, electrodynamics, etc.) in Chapter 1. Developed parallel software uses OpenMP technology (3, 9 or 27 threads) for shared memory computers. 

FORTRAN language subroutines for parallel solving 3D boundary value problems by block Gauss-Seidel method (Vanka-type smoother) are given in Chapters 2. These subroutines implement the algebraic parallelism of RMT (parallel smoothing on the finest grid).

The parallel multigrid module `RMT_3D_2020_OpenMP` is described in detail in Chapter 3. The module contains subroutines for parallel implementation of the problem-independent components of the RMT for the 3D coarse grid smoothing (geometric parallelism of RMT). 

Chapter 4 represents parallel software for solving 3D Dirichlet boundary value problem for Poisson equation (combination of geometric and algebraic parallelisms) and parallel RMT-based algorithm for solving the initial-boundary value problems.

Potential readers are graduate students and researchers working in applied and numerical mathematics as well as multigrid practitioners and software programmers for modeling physical and chemical processes in aviation and space industries, power engineering, chemical technology and other branches of mechanical engineering.

*The activity is a part of the work “Supercomputer modeling of hypervelocity impact on artificial space objects and Earth planet” supported by Russian Science Foundation (project no. 21-72-20023).*

## Mathematical background of the Robust Multigrid Technique

All the required mathematical background is provided in the book: S. I. Martynenko “The Robust Multigrid Technique: For Black-Box Software”, de Gruyter, Berlin, 2017 (https://www.degruyter.com/view/title/527481).

This book presents a detailed description of a robust pseudomultigrid algorithm for solving (initial-)boundary value problems on structured grids in a black-box manner. To overcome the problem of robustness, the presented Robust Multigrid Technique (RMT) is based on the application of the essential multigrid principle in a single grid algorithm. It results in an extremely simple, very robust and highly parallel solver with close-to-optimal algorithmic complexity and the least number of problem-dependent components. Topics covered include an introduction to the mathematical principles of multigrid methods, a detailed description of RMT, results of convergence analysis and complexity, possible expansion on unstructured grids, numerical experiments and a brief description of multigrid software, parallel RMT and estimations of speed-up and efficiency of the parallel multigrid algorithms, and finally applications of RMT for the numerical solution of the incompressible Navier-Stokes equations.

## Sequential Software for the Robust Multigrid Technique

This book (https://github.com/simartynenko/Robust_Multigrid_Technique_2020) covers sequential software for the Robust Multigrid Technique (RMT) and presents basic concepts of modern numerical methods for mathematical modeling of physical and chemical processes in the computational continuum mechanics (thermal conductivity, chemical hydrodynamics, convective heat transfer, electrodynamics, etc.). FORTRAN language subroutines for solving the boundary value problems of the computational continuum mechanics by point and block Seidel method (Vanka-type smoother) are given in Chapters 1-4. The multigrid module `RMT_3D_2020` is described in detail in Chapter 2 of this book. The module contains subroutines for implementation of the problem-independent components of the Robust Multigrid Technique for solving the three-dimensional boundary value problems. Examples of the multigrid module `RMT_3D_2020` application in static and dynamic cycles are given in Chapters 3 and 4. Chapter “Short history of RMT” represents the historical development of the Robust Multigrid Technique.

## List of modules

Multigrid modules `RMT_3D_2020_OpenMP` (chapter 3) includes subprograms for main components of RMT for parallel solving 3D boundary value problems. The multigrid module is ~~available as a static library only for Intel Fortran Compiler under Windows and GFortran under Linux operation systems~~.

Module `RMT_3D_UFG` includes the subprogram `U_Finest_Grid` for the uniform grid generation (§ 1, § 3 chapter 2 in book «Sequential Software for the Robust Multigrid Technique»).

## List of example programs

1. Program `Parallel_Finest_Grid_Smoothing_M.f90` (§ 2 chapter II): parallel solution of discrete Poisson equation (2.1) with the boundary conditions (2.4) by multicolor Vanka-type method. 

1. Program  `Integrals_M.f90` (§ 2 chapter II): parallel approximation of the triple integrals on the multigrid structure. 

1. Program  `Parallel_Coarse_Grid_Smoothing_M.f90` (§ 2 chapter II): multigrid iterations on the multigrid structures generated by the first-level dynamic grids (geometric parallelism).
