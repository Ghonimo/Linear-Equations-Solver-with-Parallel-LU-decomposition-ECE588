# Parallel LU Decomposition with POSIX Threads

## Overview
This project aims to enhance the computational efficiency of solving systems of linear equations, a cornerstone task in scientific computing, pivotal across various engineering, physics, and finance applications. By parallelizing the LU Decomposition algorithm and leveraging POSIX threads, we achieve significant speed improvements on multicore systems without compromising solution accuracy.

## Motivation
The traditional approach to solving large and complex linear systems encounters significant computational challenges, particularly concerning execution time and numerical stability. This project explores the development and implementation of a parallelized LU Decomposition algorithm to address these challenges, with a focus on leveraging the inherent potential for parallelization in modern multicore processors.

## Features
- **Parallel LU Decomposition:** Decomposes a matrix into lower (L) and upper (U) triangular matrices in a parallel fashion to simplify the solution process.
- **POSIX Thread Utilization:** Leverages POSIX threads to distribute computational tasks across multiple processors, reducing overall execution time.
- **Numerical Stability Enhancement:** Incorporates partial pivoting to increase the numerical stability of the algorithm, making it more reliable for ill-conditioned matrices.
- **Performance Optimization:** Utilizes the Bulk Synchronous Parallel (BSP) model to minimize communication overhead and ensure balanced workload distribution among processors.

## Installation
To compile and run the parallel LU Decomposition algorithm, follow these steps:

```bash
gcc -o parallel_solver parallel_solver.c -lpthread -std=gnu11
./parallel_solver
```

## Usage
After compiling, you can run the program as shown above. Modify the parallel_solver.c file to adjust the matrix size, number of threads, and other parameters as needed for your specific use case.

## Development and Debugging Tools
- Matrix Generator Script: Automates the creation of test matrices of specified dimensions.
- Data Visualization Script: Facilitates the analysis and presentation of performance data.
- Benchmark Collection Scripts: Streamlines the collection of performance benchmarks.
- Makefile for Compilation and Management: Simplifies the build process and supports clean-up operations.

## Challenges and Lessons Learned
The development of this parallel algorithm revealed several challenges, notably related to race conditions, maintaining algorithmic integrity, and ensuring numerical stability. Key lessons include the = - importance of a structured parallelization approach, the utility of debugging and optimization tools, and the critical role of numerical stability in numerical algorithms.

## Future Directions
Potential avenues for extending this work include exploring more sophisticated parallelization strategies, applying the algorithm to real-world problems, and investigating its integration into high-performance computing frameworks.