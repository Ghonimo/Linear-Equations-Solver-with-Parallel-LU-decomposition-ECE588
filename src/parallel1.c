
#include "../include/matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>

#include <math.h>

typedef struct {
    int threadId;
    int startRow;
    int endRow;
    double **A, **L, **U;
    int n;
} ThreadData;

void* parallelLUDecomposition(void* arg);

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

/////// Parallel LU Decomposition -----------------------------------
    int numThreads = 1; // Or any other way to determine the number of threads
    pthread_t threads[numThreads];
    ThreadData threadData[numThreads];

    int rowsPerThread = n / numThreads;
    for (int i = 0; i < numThreads; ++i) {
        threadData[i].threadId = i;
        threadData[i].startRow = i * rowsPerThread;
        threadData[i].endRow = (i + 1) * rowsPerThread;
        threadData[i].A = A;
        threadData[i].L = L;
        threadData[i].U = U;
        threadData[i].n = n;

        if (i == numThreads - 1) {
            threadData[i].endRow = n; // Ensure the last thread covers any remaining rows
        }

        pthread_create(&threads[i], NULL, parallelLUDecomposition, &threadData[i]);
    }

    for (int i = 0; i < numThreads; ++i) {
        pthread_join(threads[i], NULL);
    }
/////// Parallel LU Decomposition -----------------------------------

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


///////
void* parallelLUDecomposition(void* args) {
    ThreadData* data = (ThreadData*)args;
    double **A = data->A, **L = data->L, **U = data->U;
    int n = data->n;
    int startRow = data->startRow, endRow = data->endRow;
    printf("Thread %d: Starting LU decomposition from row %d to %d\n", data->threadId, startRow, endRow - 1);

    for (int i = startRow; i < endRow; i++) {
        printf("Thread %d: Computing U matrix for row %d\n", data->threadId, i);
        for (int j = i; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - sum;
            printf("U[%d][%d] = %f\n", i, j, U[i][j]);

        }
        // Ensuring L[i][i] = 1 for the diagonal of L
        L[i][i] = 1.0; // This line ensures the diagonal of L is set to 1.
        printf("L[%d][%d] set to 1 (diagonal)\n", i, i);

        printf("Thread %d: Computing L matrix for column %d\n", data->threadId, i);

        for (int j = i + 1; j < n; j++) {
            double sum = 0;
            for (int k = 0; k < i; k++) {
                sum += L[j][k] * U[k][i];
            }

                if(U[i][i] == 0) { // Check to avoid division by zero
                printf("Error: Division by zero detected at U[%d][%d]. U[%d][%d] = %f\n", i, i, i, i, U[i][i]);
                return NULL; // Early exit to avoid division by zero 
            }
            L[j][i] = (A[j][i] - sum) / U[i][i];
            printf("L[%d][%d] = %f\n", j, i, L[j][i]);

        }
    }
    printf("Thread %d: Finished LU decomposition\n", data->threadId);
    return NULL;
}
