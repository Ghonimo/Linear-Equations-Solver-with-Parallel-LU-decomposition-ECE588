#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <pthread.h>

typedef struct {
    int threadId;
    int startRow;
    int endRow;
    double **A, **L, **U;
    int n;
    pthread_mutex_t *mutexes; // Add a pointer for the mutex array
    pthread_barrier_t barrier; 
} ThreadData;

void* parallelLUDecomposition(void* arg);

struct timespec StartTime;
struct timespec EndTime;

// Function to perform LU decomposition (sequential version for reference)
void luDecomposition(double** A, double** L, double** U, int n) {
    for (int i = 0; i < n; i++) {
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);
            U[i][k] = A[i][k] - sum;
        }

        for (int k = i; k < n; k++) {
            if (i == k)
                L[i][i] = 1; 
            else {
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);
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

    fscanf(file, "%d", n); 

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

void* parallelLUDecomposition(void* args) {
    ThreadData* data = (ThreadData*)args;
    double **A = data->A, **L = data->L, **U = data->U;
    int n = data->n;
    int startRow = data->startRow, endRow = data->endRow;
    //pthread_mutex_t *mutexes = data->mutexes; 
    pthread_barrier_t *barrier = &data->barrier;

    for (int i = startRow; i < endRow; i++) {
        printf("Thread %d: Starting to work on row %d of U\n", data->threadId, i);

        // Precalculate sums
        double rowSums[n];
        for (int k = i; k < n; k++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * U[j][k];
            }
            rowSums[k] = sum;
        }

        // Calculate U elements in parallel
        for (int k = i; k < n; k++) {
            double sum = 0; 
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * U[j][k];
            }
            pthread_mutex_lock(&data->mutexes[i]); // Access mutexes through data->mutexes
            printf("Thread %d: Acquired mutex %d (U update)\n", data->threadId, i);
            printf("Thread %d: Calculating U[%d][%d], relevant values of L and U: ... \n", data->threadId, i, k); // Add debugging output here
            U[i][k] = A[i][k] - sum;
            pthread_mutex_unlock(&data->mutexes[i]); // Access mutexes through data->mutexes
            printf("Thread %d: Released mutex %d (U update)\n", data->threadId, i);
        }

        printf("Thread %d: Finished updating row %d of U\n", data->threadId, i);
        pthread_barrier_wait(barrier);
        printf("Thread %d: Passed the barrier\n", data->threadId, i);

        // Calculate L elements 
        for (int k = i + 1; k < n; k++) { 
            double sum = 0; 
            for (int j = 0; j < i; j++) { 
                sum += L[k][j] * U[j][i];
            } 
            pthread_mutex_lock(&data->mutexes[k]); // Use data->mutexes
            printf("Thread %d: Acquired mutex %d (L update)\n", data->threadId, k);
            L[k][i] = (A[k][i] - sum) / U[i][i];
            pthread_mutex_unlock(&data->mutexes[k]); // Use data->mutexes
            printf("Thread %d: Released mutex %d (L update)\n", data->threadId, k);
        }

        printf("Thread %d: Finished updating row %d of L\n", data->threadId, i);
    }
    return NULL;
}


// Main Function
int main(int argc, char* argv[]) {
    if (argc < 2 || argc > 3) {
        printf("Usage: %s <matrix_file> [output_file]\n", argv[0]);
        return EXIT_FAILURE;
    }

    int ret;
    int n; 
    double **A, **L, **U, *b, *y, *x;
    struct timespec start, end; 

    clock_gettime(CLOCK_REALTIME, &start); 

    readMatrixFromFile(argv[1], &A, &b, &n);

    L = (double**)malloc(n * sizeof(double*));
    U = (double**)malloc(n * sizeof(double*));
    y = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        L[i] = (double*)malloc(n * sizeof(double));
        U[i] = (double*)malloc(n * sizeof(double));
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            L[i][j] = 0;
            U[i][j] = 0;
        }
    }

    ret = clock_gettime(CLOCK_REALTIME, &StartTime);
    assert(ret == 0);

    int numThreads = 4; 
    pthread_t threads[numThreads];
    ThreadData threadData[numThreads];

    pthread_mutex_t **threadMutexes = malloc(numThreads * sizeof(pthread_mutex_t*));
    for (int i = 0; i < numThreads; i++) {
        threadMutexes[i] = malloc(n * sizeof(pthread_mutex_t));
        for (int j = 0; j < n; j++) {
            pthread_mutex_init(&threadMutexes[i][j], NULL);
        }
    }

    pthread_barrier_init(&threadData[0].barrier, NULL, numThreads);

    int rowsPerThread = n / numThreads;
    for (int i = 0; i < numThreads; i++) {
        threadData[i].threadId = i;
        threadData[i].startRow = i * rowsPerThread;
        threadData[i].endRow = (i + 1) * rowsPerThread;
        threadData[i].A = A;
        threadData[i].L = L;
        threadData[i].U = U;
        threadData[i].n = n;
        threadData[i].mutexes = threadMutexes[i]; 

        if (i == numThreads - 1) {
            threadData[i].endRow = n; 
        }

        pthread_create(&threads[i], NULL, parallelLUDecomposition, &threadData[i]);
    }

    for (int i = 0; i < numThreads; i++) { 
        pthread_join(threads[i], NULL);
    }

    for (int i = 0; i < numThreads; i++) {
        for (int j = 0; j < n; j++) {
            threadData[i].mutexes = threadMutexes[i]; // Assign the mutex array 
        }
        free(threadMutexes[i]);
    }
    free(threadMutexes); 

    forwardSubstitution(L, b, y, n);
    backwardSubstitution(U, y, x, n);

    ret = clock_gettime(CLOCK_REALTIME, &EndTime);
    assert(ret == 0);

    clock_gettime(CLOCK_REALTIME, &end); 
    double time_taken = end.tv_sec - start.tv_sec + (end.tv_nsec - start.tv_nsec) / 1e9;

    FILE* outputFile;
    if (argc == 3) {
        outputFile = fopen(argv[2], "w");
        if (outputFile == NULL) {
            perror("Error opening output file");
            freeMemory(A, L, U, b, y, x, n);
            return EXIT_FAILURE;
        }
    } else {
        outputFile = stdout; 
    }

    fprintf(outputFile, "Solution: \n");
    for (int i = 0; i < n; i++) {
        fprintf(outputFile, "x[%d] = %f\n", i, x[i]);
    }
    fprintf(outputFile, "\nTime taken: %.9f seconds\n", time_taken);

    if (outputFile != stdout) {
        fclose(outputFile);
    }

    freeMemory(A, L, U, b, y, x, n);

    return 0;
}
