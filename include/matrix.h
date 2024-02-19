// matrix_.h
#ifndef MATRIX_H
#define MATRIX_H

#define N 3 // If N is fixed; otherwise, remove it from here and define it elsewhere

void luDecomposition(double A[N][N], double L[N][N], double U[N][N]);
void forwardSubstitution(double L[N][N], double b[N], double y[N]);
void backwardSubstitution(double U[N][N], double y[N], double x[N]);

#endif // MATRIX_H
