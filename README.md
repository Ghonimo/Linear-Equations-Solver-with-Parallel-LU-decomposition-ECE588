# Parallel Matrix Solver 🚀

Welcome to the Parallel Matrix Solver project! This project is designed to solve large systems of linear equations using LU Decomposition, leveraging the power of parallel computing to speed up computations. Perfect for educational purposes and anyone interested in exploring the efficiency of parallel algorithms.

## Project Directory Structure 📚

Dive into the heart of our Parallel Matrix Solver project. Here's how everything is organized:

### Source Code 🧬
`src/` - This is where the magic happens! Our source code directory includes both serial and parallel implementations of the LU Decomposition algorithm.
- `parallel_solver.c`: Embark on a journey with our parallel version, optimized with Pthreads for concurrency.
- `sequential_solver.c`: The classic approach to LU Decomposition. Perfect for understanding the basics and performance comparison.
- `...`: Other variations and utilities to explore different implementations of LU Decomposition.

### Headers 📑
`include/` - All our project-wide declarations live here.
- `matrix.h`: Essential definitions and function prototypes for matrix operations.

### Documentation 📖
`docs/` - Need guidance? Our docs have got you covered!
- `setup.md`: Get up and running with detailed setup instructions.
- `usage.md`: Learn how to wield our solver with comprehensive usage instructions.
- `...`: Additional documentation and design logs for our project.

### Benchmarks 📈
`benchmarks/` - Curious about performance? Check out our benchmarks!
- `results.md`: Dive into detailed performance analysis and see how our parallel implementation stacks up against the serial one.

### Matrices 🧩
`matrices/` - Sample matrices for testing and getting a feel of the solver's prowess.
- `...`: A variety of matrices to challenge and benchmark our solver.

### Solutions 🏁
`matrices_solution/` - Wondering if you got it right? Here are the solutions for our sample matrices.
- `...`: Solutions for all provided sample matrices to verify your results.

### Utilities 🔧
`utilities/` - Our utility belt for benchmarking and matrix generation.
- `...`: Tools to help you generate matrices, run benchmarks, and more.

```
## Prerequisites

* A C compiler (GCC, Clang, etc.)
* Pthreads library 

## Setup

1. Clone the repository:
```bash
   git clone git@github.com:Ghonim/ECE588_Parallel_Matrix.git
```

2. Navigate to the project directory:
```bash
cd ECE588_Parallel_Matrix
```

[Setup Instructions](docs/setup.md)

[Usage Instructions](docs/usage.md)

## Building the Project

The project includes a Makefile for easy compilation of the source code into executable binaries. You can use the following commands to build the project:

- To build all versions of the solver (sequential, parallel, and parallel with pivoting):
  ```bash
  make all
   ```
- To build specific versions of the solver, use one of the following commands:
For the sequential version:
   ```bash
   make sequential
   ```
- For the parallel version without pivoting:
   ```bash
   make parallel
   ```
- For the parallel version with pivoting:
   ```bash
   make pivoting
   ```

