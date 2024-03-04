// main.c
// Author: Mohamed Ghonim
// Created: 02/18/2024
// Last Modified: 02/27/2024

// Debug: Need to fix the parallel version of the code, it runs but the output is not correct. 

#include "../include/matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>

#define MAX_THREADS 32

struct timespec StartTime;
struct timespec EndTime;

typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t cond;
    int count;
    int num_threads;
} my_barrier_t;

int my_barrier_init(my_barrier_t *barrier, unsigned int num_threads) {
    pthread_mutex_init(&barrier->mutex, NULL);
    pthread_cond_init(&barrier->cond, NULL);
    barrier->count = 0;
    barrier->num_threads = num_threads;
    return 0; 
}


int my_barrier_destroy(my_barrier_t *barrier) {
    pthread_mutex_destroy(&barrier->mutex);
    pthread_cond_destroy(&barrier->cond);
    return 0;
}

int my_barrier_wait(my_barrier_t *barrier) {
    pthread_mutex_lock(&barrier->mutex);
    barrier->count++;
    if (barrier->count == barrier->num_threads) {
        barrier->count = 0; // Reset for next use
        pthread_cond_broadcast(&barrier->cond);
    } else {
        while (pthread_cond_wait(&barrier->cond, &barrier->mutex) != 0);
    }
    pthread_mutex_unlock(&barrier->mutex);
    return 0;
}

int n;                      // Matrix size
double **A, **L, **U, *b, *y, *x;
int numThreads;             // Number of threads for parallel computation
my_barrier_t barrier;  // For synchronizing threads after LU decomposition

// Thread function for LU decomposition
void* luDecompositionThread(void* arg) {
    long threadId = (long)arg;
    int startRow = threadId * n / numThreads;
    int endRow = (threadId + 1) * n / numThreads;

    for (int i = startRow; i < endRow; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * U[j][k];
            }
            U[i][k] = A[i][k] - sum;
        }
        for (int k = i + 1; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[k][j] * U[j][i];
            }
            L[k][i] = (A[k][i] - sum) / U[i][i];
        }
    }

    my_barrier_wait(&barrier); // Synchronize threads after LU decomposition

    return NULL;
}

// Main Function
int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        printf("Usage: %s <matrix_file> [output_file]\n", argv[0]);
        return EXIT_FAILURE;
    }

    int ret;
    struct timespec start, end; // Timing variables

    clock_gettime(CLOCK_REALTIME, &start); // Record start time

    // Read matrix A and vector b from file
    readMatrixFromFile(argv[1], &A, &b, &n);

    // Parse number of threads from command-line arguments
    int numThreads = (argc == 3) ? atoi(argv[2]) : 4; // Default to 4 threads if not specified
    if (numThreads <= 0 || numThreads > MAX_THREADS) numThreads = MAX_THREADS;

    // Allocate memory for L, U, y, and x
    L = (double**)malloc(n * sizeof(double*));
    U = (double**)malloc(n * sizeof(double*));
    y = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        L[i] = (double*)malloc(n * sizeof(double));
        U[i] = (double*)malloc(n * sizeof(double));
    // Initialize L and U to zeros
        for (int j = 0; j < n; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    // Initialize barrier
    my_barrier_init(&barrier, numThreads);

    // Starting the clock
    struct timespec StartTime, EndTime;
    ret = clock_gettime(CLOCK_REALTIME, &StartTime);
    assert(ret == 0);

    // Create threads for LU decomposition
    pthread_t threads[MAX_THREADS];
    for (long i = 0; i < numThreads; i++) {
        pthread_create(&threads[i], NULL, luDecompositionThread, (void*)i);
    }

    // Join threads after LU decomposition
    for (int i = 0; i < numThreads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Destroy barrier
    my_barrier_destroy(&barrier);

    // Ending Clock
    ret = clock_gettime(CLOCK_REALTIME, &EndTime);
    assert(ret == 0);
    clock_gettime(CLOCK_REALTIME, &end); // Record end time

    // Perform forward and backward substitution sequentially
    forwardSubstitution(L, b, y, n);
    backwardSubstitution(U, y, x, n);

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


///////////

