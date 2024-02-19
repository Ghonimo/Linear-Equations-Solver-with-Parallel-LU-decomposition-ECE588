# Parallel Equation Solver

This project implements a parallel linear equation solver for large matrices, employing both serial and parallel LU Decomposition algorithms using Pthreads. It explores performance enhancement achieved through multicore execution strategies.

## Project Structure
```bash
Parallel Equation Solver
│
├── src/                    # Source files
│   ├── main.c              # Main program entry point
│   ├── lu_serial.c         # Serial LU Decomposition implementation
│   ├── lu_serial.h         # Header for serial LU Decomposition
│   ├── lu_parallel.c       # Parallel LU Decomposition implementation using Pthreads
│   ├── lu_parallel.h       # Header for parallel LU Decomposition
│   └── utilities.c         # Utility functions for matrix operations and timing
│
├── include/                # Header files
│   ├── matrix.h            # Definitions and functions for matrix operations
│   └── timing.h            # Timing and performance measurement utilities
│
├── docs/                   # Documentation files
│   ├── setup.md            # Setup instructions
│   ├── usage.md            # Usage instructions
│   └── development.md      # Notes on development decisions and project structure
│
├── tests/                  # Test files
│   ├── test_serial.c       # Tests for serial implementation
│   └── test_parallel.c     # Tests for parallel implementation
│
├── benchmarks/             # Benchmark scripts and results
│   ├── benchmark_script.sh # Script to run benchmarks
│   └── results.md          # Benchmark results and analysis
│
├── Makefile                # Makefile for building the project
└── README.md               # Project overview and general instructions
```
## Prerequisites

* A C compiler (GCC, Clang, etc.)
* Pthreads library 

## Setup

1. Clone the repository:
```bash
   git clone https://github.com/Ghonimo/ECE588_Parallel_Matrix
```

2. Navigate to the project directory:
```bash
cd ECE588_Parallel_Matrix
```

[Setup Instructions](docs/setup.md)

[Usage Instructions](docs/usage.md)
