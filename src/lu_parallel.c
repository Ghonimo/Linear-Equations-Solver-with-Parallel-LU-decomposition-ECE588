// lu_parallel.c
// Author: Phanindra Vemireddy
// Created: 02/23/2024
// Last Modified: 02/23/2024


#include "../include/matrix.h"
#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <sys/time.h>
#include <pthread.h>
#include <errno.h>
#include <sys/errno.h>
#include <time.h>

#define M_NPROCS 16
int NumThreads;
int Count;
int Remainder;
struct timespec StartTime;
struct timespec EndTime;

pthread_mutex_t SyncLock;
pthread_cond_t SyncCV;
int SyncCount;

double **A, **L, **U, *b, *y, *x;
int n; // Matrix size

/*Barrier function is same from synchronization from sum.c*/

void Barrier()
{
 int ret;
 pthread_mutex_lock(&SyncLock);
 SyncCount++;
 if(SyncCount == NumThreads)
 {
 ret = pthread_cond_broadcast(&SyncCV);
 assert(ret == 0);
 SyncCount=0;
 }
 else
 {
 ret = pthread_cond_wait(&SyncCV, &SyncLock);
 assert(ret == 0);
 }
 pthread_mutex_unlock(&SyncLock);
}

// Function to perform LU decomposition
// void luDecomposition_p(double** A, double** L, double** U, int start_row, int end_row, int start_col, int end_col) {
void luDecomposition_p(int start_row, int end_row, int start_col, int end_col) {
        printf("DEBUG: luDecomposition_p entered with start_row=%d, end_row=%d, start_col=%d, end_col=%d\n", start_row, end_row, start_col, end_col);

    for (int i = start_row; i < end_row; i++) {
        for (int k = i; k < end_col; k++) {
            // Summation of L(i, j) * U(j, k)
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (L[i][j] * U[j][k]);

            // Evaluating U(i, k)
            U[i][k] = A[i][k] - sum;
            printf("DEBUG: Inside outer loop, i=%d, k=%d\n", i, k);
        }

        for (int k = i; k < end_col; k++) {
            if (i == k)
                L[i][i] = 1; // Diagonal as 1
            else {
                // Summation of L(k, j) * U(j, i)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (L[k][j] * U[j][i]);

                // Evaluating L(k, i)
                L[k][i] = (A[k][i] - sum) / U[i][i];
                printf("DEBUG: Inside inner loop, i=%d, k=%d\n", i, k);

            }
        }
    }
    printf("DEBUG: luDecomposition_p exiting\n"); 

}

void* thread_task(void* arg) {
    long int thread_id = (long int) arg;
    int block_size = n / NumThreads;
    int start_row = (thread_id / 2) * block_size;
    int start_col = (thread_id % 2) * block_size;
    int end_row = start_row + block_size;
    int end_col = start_col + block_size;

    int s_row, e_row, s_col, e_col;
    for (int i = start_row; i < n; i += 2 * block_size) {
        for (int j = start_col; j < n; j += 2 * block_size) {
            // for (int x = i; x < i + block_size && x < n; x++) {
            //     for (int y = j; y < j + block_size && y < n; y++) {
            //         result_matrix[x][y] = *thread_id;
            //     }
            // }
            s_row = i;
            e_row = i + block_size;
            s_col = j;
            e_col = j + block_size;
            // luDecomposition_p(A, L, U, s_row, e_row, s_col, e_col);
            luDecomposition_p(s_row, e_row, s_col, e_col);
        }
    }

    pthread_exit(NULL);
}

void initializeLU(double*** L, double*** U, int n) {
    *L = (double**)malloc(n * sizeof(double*));
    *U = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        (*L)[i] = (double*)malloc(n * sizeof(double));
        (*U)[i] = (double*)malloc(n * sizeof(double));
    }
}


// Main Function
int main(int argc, char* argv[]) {
    // if (argc < 2 || argc > 5) {
    //     printf("Invalid Usage: %s <matrix_file> [output_file]\n", argv[0]);
    //     return EXIT_FAILURE;
    // }

    int ret;
    struct timespec start, end; // Timing variables

    /*Loop Variables*/
    pthread_t* Threads;
    pthread_attr_t attr;
    long int ThreadID;

    
    /*Arguments passed by command line*/

    // if(argc<2)
    // {
    // fprintf(stderr, "Syntax: %s <numProcesors>\nExiting Program...\n", argv[0]);
    // exit(1);
    // }
    NumThreads = 1;
    if (NumThreads < 1 || NumThreads >M_NPROCS)
    {
    fprintf(stderr,"Number of processors has to be between 1 and %d\nExiting Program...\n",M_NPROCS);
    exit(1);
    }

    /*Initializing threads*/

    Threads = (pthread_t *) malloc(sizeof(pthread_t) * NumThreads);
    assert(Threads!=NULL);

    /*Step Calculation*/

    // Count = (Xmax-1) / NumThreads;
    // Remainder = (Xmax-1) % NumThreads;

    /* Initialize thread attribute */

    pthread_attr_init(&attr);
    pthread_attr_setscope(&attr, PTHREAD_SCOPE_SYSTEM);
    ret=pthread_mutex_init(&SyncLock,NULL);
    assert(ret==0);
    ret=pthread_cond_init(&SyncCV,NULL);
    assert(ret==0);
    SyncCount=0;

    /*Starting the clock*/
    ret = clock_gettime(CLOCK_REALTIME, &StartTime);
    assert(ret == 0);


    
    y = (double*)malloc(n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));

    initializeLU(&L, &U, n);
    
    
    // Read matrix A and vector b from file
    // readMatrixFromFile(argv[2], &A, &b, &n);
    readMatrixFromFile("matrices/200x200_1.txt", &A, &b, &n);

    // Starting the clock
    ret = clock_gettime(CLOCK_REALTIME, &StartTime);
    assert(ret == 0);

     /*Start threading*/

    for(ThreadID=0; ThreadID < NumThreads; ThreadID++)
    {
    ret = pthread_create(&Threads[ThreadID], &attr, thread_task, (void*) ThreadID);
    assert(ret == 0);
    }

    /*Join threads after termination*/

    for(ThreadID=0; ThreadID < NumThreads; ThreadID++)
    {
    ret = pthread_join(Threads[ThreadID], NULL);
    assert(ret == 0);
    }

    /*Ending Clock*/

    ret = clock_gettime(CLOCK_REALTIME, &EndTime);
    assert(ret == 0);


    // Perform LU decomposition and solve the system
    forwardSubstitution(L, b, y, n);
    backwardSubstitution(U, y, x, n);

    // // Print solution
    // printf("Solution: \n");
    // for (int i = 0; i < n; i++) {
    //     printf("x[%d] = %f\n", i, x[i]);
    // }


    // Write solution to file if specified
    FILE* outputFile;
    if (argc == 4) {
        outputFile = fopen(argv[3], "w");
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
    //printf("\nTime = %lld nanoseconds\t(%ld.%09ld sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);
    printf("\nTime = %llu nanoseconds\t(%llu.%09llu sec)\n", runtime, runtime / 1000000000, runtime % 1000000000);

 
    fprintf(outputFile, "Solution: \n");
    for (int i = 0; i < n; i++) {
        fprintf(outputFile, "x[%d] = %f\n", i, x[i]);
    }
    // fprintf(outputFile, "\nTime taken: %.9f seconds\n", time_taken);
    
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
