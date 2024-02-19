// matrix_.h
#ifndef MATRIX_H
#define MATRIX_H

#define N 3 // If N is fixed; otherwise, remove it from here and define it elsewhere

// Function Prototypes
void luDecomposition(double** A, double** L, double** U, int n);
void forwardSubstitution(double** L, double* b, double* y, int n);
void backwardSubstitution(double** U, double* y, double* x, int n);
void readMatrixFromFile(const char* filename, double*** A, double** b, int* n);
void freeMemory(double** A, double** L, double** U, double* b, double* y, double* x, int n);

#endif // MATRIX_H
