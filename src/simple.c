// Author: Mohamed Ghonim
// Created: 02/18/2024
// Last Modified: 02/18/2024
// Functionality: Perform LU decomposition and solve a system of linear equations using forward and backward substitution.
// Version: 0.2 : readin the matrix from a file passed as an argument to the program

    // double A[N][N] = {{1, 1, -1}, {1, -2, 3}, {2, 3, 1}};
    // double b[N] = {4, -6, 7};

#include "../include/matrix.h"
#include <stdio.h>
#include <stdlib.h> // For exit() and EXIT_FAILURE

#define N 3 // Size of the matrix (3x3)

// Assume LU Decomposition, Forward Substitution, and Backward Substitution
void readMatrixFromFile(const char* filename, double A[N][N], double b[N]) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (fscanf(file, "%lf", &A[i][j]) != 1) {
                perror("Error reading matrix from file");
                fclose(file);
                exit(EXIT_FAILURE);
            }
        }
    }

    for (int i = 0; i < N; i++) {
        if (fscanf(file, "%lf", &b[i]) != 1) {
            perror("Error reading vector from file");
            fclose(file);
            exit(EXIT_FAILURE);
        }
    }

    fclose(file);
}

///////
// Function to perform LU decomposition
void luDecomposition(double A[N][N], double L[N][N], double U[N][N]) {
    int i, j, k;
    for (i = 0; i < N; i++) {
        for (j = 0; j < N; j++) {
            if (j < i)
                L[j][i] = 0;
            else {
                L[j][i] = A[j][i];
                for (k = 0; k < i; k++) {
                    L[j][i] = L[j][i] - L[j][k] * U[k][i];
                }
            }
        }
        for (j = 0; j < N; j++) {
            if (j < i)
                U[i][j] = 0;
            else if (j == i)
                U[i][j] = 1;
            else {
                U[i][j] = A[i][j] / L[i][i];
                for (k = 0; k < i; k++) {
                    U[i][j] = U[i][j] - ((L[i][k] * U[k][j]) / L[i][i]);
                }
            }
        }
    }
}

// Function to solve the equation Ly = b
void forwardSubstitution(double L[N][N], double b[N], double y[N]) {
    for (int i = 0; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= L[i][j] * y[j];
        }
        y[i] = y[i] / L[i][i];
    }
}

// Function to solve the equation Ux = y
void backwardSubstitution(double U[N][N], double y[N], double x[N]) {
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= U[i][j] * x[j];
        }
        // No division by U[i][i] since U[i][i] = 1
    }
}

//

int main(int argc, char* argv[]) {
    if (argc != 2) {
        printf("Usage: %s <matrix_file>\n", argv[0]);
        return EXIT_FAILURE;
    }

    double A[N][N];
    double b[N];
    double L[N][N] = {0};
    double U[N][N] = {0};
    double y[N] = {0};
    double x[N] = {0};

    // Read matrix A and vector b from file
    readMatrixFromFile(argv[1], A, b);

    luDecomposition(A, L, U);
    forwardSubstitution(L, b, y);
    backwardSubstitution(U, y, x);

    printf("Solution: \n");
    for (int i = 0; i < N; i++) {
        printf("x[%d] = %f\n", i, x[i]);
    }

    return 0;
}
