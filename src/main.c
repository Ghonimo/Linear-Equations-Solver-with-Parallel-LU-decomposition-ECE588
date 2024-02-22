// main.c
// Author: Mohamed Ghonim
// Created: 02/18/2024
// Last Modified: 02/18/2024
// Functionality: Perform LU decomposition and solve a system of linear equations using forward and backward substitution.
// Version: 0.3 : read in the matrix from a file passed as an argument to the program
//              : added a function to free the allocated memory
//              : the matrix size is parameterized, and passed as the first line in the matrix file

// Version: 0.4 : allow writing the solution to a file passed as an argument to the program
//              : if no file is passed, the solution will be printed in the terminal
//              : example usage: matrices/py_generated/5000x5000.txt matrices_solution/5000x5000.txt


#include "../include/matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

struct timespec StartTime;
struct timespec EndTime;

// Main Function
int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        printf("Usage: %s <matrix_file> [output_file]\n", argv[0]);
        return EXIT_FAILURE;
    }

    int ret;
    int n; // Matrix size
    double **A, **L, **U, *b, *y, *x;
    struct timespec start, end; // Timing variables

    clock_gettime(CLOCK_REALTIME, &start); // Record start time

    // Read matrix A and vector b from file
    readMatrixFromFile(argv[1], &A, &b, &n);

    // Allocate memory for L, U, y, and x
    L = (double**)malloc(n * sizeof(double*));
    U = (double**)malloc(n * sizeof(double*));
    y = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        L[i] = (double*)malloc(n * sizeof(double));
        U[i] = (double*)malloc(n * sizeof(double));
    }

    // Initialize L and U to zeros
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    // Starting the clock
    struct timespec StartTime, EndTime;
    ret = clock_gettime(CLOCK_REALTIME, &StartTime);
    assert(ret == 0);

    // Perform LU decomposition and solve the system
    luDecomposition(A, L, U, n);
    forwardSubstitution(L, b, y, n);
    backwardSubstitution(U, y, x, n);

    // Ending Clock
    ret = clock_gettime(CLOCK_REALTIME, &EndTime);
    assert(ret == 0);
    // // Print solution
    // printf("Solution: \n");
    // for (int i = 0; i < n; i++) {
    //     printf("x[%d] = %f\n", i, x[i]);
    // }

    clock_gettime(CLOCK_REALTIME, &end); // Record end time
    double time_taken = end.tv_sec - start.tv_sec + (end.tv_nsec - start.tv_nsec) / 1e9; // Calculate elapsed time in seconds

    // Write solution to file if specified
    FILE* outputFile;
    if (argc == 3) {
        outputFile = fopen(argv[2], "w");
        if (outputFile == NULL) {
            perror("Error opening output file");
            freeMemory(A, L, U, b, y, x, n);
            return EXIT_FAILURE;
        }
    } else {
        outputFile = stdout; // Use standard output if no file is specified
    }

    // Use outputFile instead of printf for output
    // fprintf(outputFile, "Solution: \n");
    // for (int i = 0; i < n; i++) {
    //     fprintf(outputFile, "x[%d] = %f\n", i, x[i]);
    // }
    
    unsigned long long int runtime = 1000000000 * (EndTime.tv_sec - StartTime.tv_sec) + EndTime.tv_nsec - StartTime.tv_nsec;
    printf("\nTime = %lld nanoseconds\t(%ld.%09ld sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
 
    fprintf(outputFile, "Solution: \n");
    for (int i = 0; i < n; i++) {
        fprintf(outputFile, "x[%d] = %f\n", i, x[i]);
    }
    fprintf(outputFile, "\nTime taken: %.9f seconds\n", time_taken);
    
    if (outputFile != stdout) {
        fclose(outputFile);
    }

    // Free allocated memory
    freeMemory(A, L, U, b, y, x, n);

    return 0;
}

// Function to perform LU decomposition
void luDecomposition(double** A, double** L, double** U, int n) {
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            // Summation of L(i, j) * U(j, k)
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);

            // Evaluating U(i, k)
            U[i][k] = A[i][k] - sum;
        }

        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal as 1
            else {
                // Summation of L(k, j) * U(j, i)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);

                // Evaluating L(k, i)
                L[k][i] = (A[k][i] - sum) / U[i][i];
            }
        }
    }
}

// Function to solve the equation Ly = b for y
void forwardSubstitution(double** L, double* b, double* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int k = 0; k < i; k++)
            y[i] -= L[i][k] * y[k];
        y[i] = y[i] / L[i][i];
    }
}

// Function to solve the equation Ux = y for x
void backwardSubstitution(double** U, double* y, double* x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}


// Function to read matrix from file
void readMatrixFromFile(const char* filename, double*** A, double** b, int* n) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d", n); // Read the size of the matrix from the first line

    *A = (double**)malloc(*n * sizeof(double*));
    *b = (double*)malloc(*n * sizeof(double));
    for (int i = 0; i < *n; i++) {
        (*A)[i] = (double*)malloc(*n * sizeof(double));
        for (int j = 0; j < *n; j++) {
            fscanf(file, "%lf", &(*A)[i][j]);
        }
    }

    for (int i = 0; i < *n; i++) {
        fscanf(file, "%lf", &(*b)[i]);
    }

    fclose(file);
}

// Function to free dynamically allocated memory
void freeMemory(double** A, double** L, double** U, double* b, double* y, double* x, int n) {
    for (int i = 0; i < n; i++) {
        free(A[i]);
        free(L[i]);
        free(U[i]);
    }
    free(A);
    free(L);
    free(U);
    free(b);
    free(y);
    free(x);
}
