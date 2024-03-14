// cc -lpthread src/parallel_solver.c -o bin/parallel_solver -std=gnu11

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

double **a, **a_duplicate, **u, **l, **permutation_matrix, **PA, **LU;
double *b, *x, *y;  // b is the last line of the matrix file, x and y are solution vectors
int n;
int numThreads;     // Number of threads to be used

struct timespec start, end; // Timing variables
struct thread_data {
    int id;         // thread ID
    int k;          // current step of the decomposition
};

/**
 * Function to solve the equation Ly = b for y
 *  L The lower triangular matrix
 *  b The right-hand side vector
 *  y The solution vector
 * n The size of the matrix
 */
void forwardSubstitution(double** L, double* b, double* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int k = 0; k < i; k++)
            y[i] -= L[i][k] * y[k];
        y[i] = y[i] / L[i][i];
    }
}

/**
 * Function to solve the equation Ux = y for x
 * U The upper triangular matrix
 * y The right-hand side vector
 * x The solution vector
 * n The size of the matrix
 */
void backwardSubstitution(double** U, double* y, double* x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

/**
 * Function to print a matrix
 * matrix The matrix to be printed
 * n The size of the matrix
 */
void print(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f  ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

/**
 * Function to multiply two matrices
 * matrix1 The first matrix
 * matrix2 The second matrix
 * product The result matrix
 * n The size of the matrices
 */
void matrix_multiply(double** matrix1, double** matrix2, double** product, int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            product[i][j] = 0;
            for (int k = 0; k < n; ++k) {
                product[i][j] += matrix1[i][k] * matrix2[k][j];
            }
        }
    }
}

/**
 * Function to be executed by each thread to assign values to matrix A in parallel
 * thread_data The thread data containing the thread ID and the value of k
 */
void* parallel_portion(void* thread_data) {
    struct thread_data* my_data;
    my_data = (struct thread_data*) thread_data;
    int id = my_data->id;
    int k = my_data->k;

    printf("Thread %d: Computing step k: %d\n", id, k);
    int interation_per_thread = n - 1 - k;
    int start = (k + 1) + id * interation_per_thread / numThreads;
    int end = (k + 1) + (id + 1) * interation_per_thread / numThreads < n ? (k + 1) + (id + 1) * interation_per_thread / numThreads : n;
    for (int i = start; i < end; i++) {
        for (int j = k + 1; j < n; j++) {
            a[i][j] -= l[i][k] * u[k][j];
            printf("Thread %d: A[%d][%d] = %f\n", id, i, j, a[i][j]);
        }
    }
    pthread_exit(NULL); 
}

/**
 * Function to read a matrix from a file
 * filename The name of the file containing the matrix
 */
void readMatrixFromFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d", &n); // Read the size of the matrix from the first line

    a = (double**)malloc(n * sizeof(double*));
    a_duplicate = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        a[i] = (double*)malloc(n * sizeof(double));
        a_duplicate[i] = (double*)malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) {
            fscanf(file, "%lf", &a[i][j]);
            a_duplicate[i][j] = a[i][j]; // Copy the value to a_duplicate
        }
    }
    
    // Assuming the vector b is in the last line after the matrix
    b = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        fscanf(file, "%lf", &b[i]);
    }

    fclose(file);
}

void performPartialPivoting(double** a, double** l, int* p, int k) {
    double max_value = 0.0;
    int max_index = k;
    for (int i = k; i < n; i++) {
        double absolute_value = fabs(a[i][k]);
        if (absolute_value > max_value) {
            max_value = absolute_value;
            max_index = i;
        }
    }

    if (max_value == 0.0) {
        fprintf(stderr, "This is a singular matrix, LU Decomposition is not possible\n");
        exit(EXIT_FAILURE);
    }

    // Swap p[k] and p[max_index]
    int temp_p = p[k];
    p[k] = p[max_index];
    p[max_index] = temp_p;

    // Swap rows in a
    double* temp_a = a[k];
    a[k] = a[max_index];
    a[max_index] = temp_a;

    // Swap rows in l for the elements before diagonal (k)
    for (int i = 0; i < k; i++) {
        double temp_l = l[k][i];
        l[k][i] = l[max_index][i];
        l[max_index][i] = temp_l;
    }
}

void freeUpMemory(double*** a, double*** a_duplicate, double*** u, double*** l, double*** permutation_matrix, double*** PA, double*** LU, int** p) {
    for (int i = 0; i < n; i++) {
        free((*a)[i]);
        free((*a_duplicate)[i]);
        free((*u)[i]);
        free((*l)[i]);
        free((*permutation_matrix)[i]);
        free((*PA)[i]);
        free((*LU)[i]);
    }
    free(*a);
    free(*a_duplicate);
    free(*u);
    free(*l);
    free(*permutation_matrix);
    free(*PA);
    free(*LU);
    free(*p);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: %s <file path> <number of numThreads>\n", argv[0]);
        return 1;
    }
    clock_gettime(CLOCK_REALTIME, &start); // Record start time

    // Read matrix from file
    readMatrixFromFile(argv[1]);
    numThreads = atoi(argv[2]);

    // Memory allocation for L, U, permutation_matrix, PA, LU
    u = (double**)malloc(n * sizeof(double*));
    l = (double**)malloc(n * sizeof(double*));
    permutation_matrix = (double**)malloc(n * sizeof(double*));
    PA = (double**)malloc(n * sizeof(double*));
    LU = (double**)malloc(n * sizeof(double*));

    x = (double*)malloc(n * sizeof(double));
    y = (double*)malloc(n * sizeof(double));

    for (int i = 0; i < n; i++) {
        u[i] = (double*)malloc(n * sizeof(double));
        l[i] = (double*)malloc(n * sizeof(double));
        permutation_matrix[i] = (double*)malloc(n * sizeof(double));
        PA[i] = (double*)malloc(n * sizeof(double));
        LU[i] = (double*)malloc(n * sizeof(double));
    }

    int* p = (int*)malloc(n * sizeof(int));
    for (int i = 0; i < n; i++) {
        p[i] = i;
        for (int j = 0; j < n; j++) { 
            if (j < i) { 
                u[i][j] = 0.0;
                l[i][j] = 0.0;
            } else if (j == i) { 
                l[i][j] = 1.0;
            } else {
                u[i][j] = 0.0;
                l[i][j] = 0.0;
            }
        }
    }

    double max_value;
    int max_index;
    pthread_t threads[numThreads];
    struct thread_data thread_data_array[numThreads];

    printf("Starting LU Decomposition with %d threads for %dx%d matrix.\n", numThreads, n, n);

    for (int k = 0; k < n; k++) {
        // Call the partial pivoting function
        //performPartialPivoting(a, l, p, k);

        u[k][k] = a[k][k];
        for (int i = k + 1; i < n; i++) {
            l[i][k] = a[i][k] / u[k][k];
            u[k][i] = a[k][i];
        }

        for (int i = 0; i < numThreads; i++) {
            thread_data_array[i].id = i;
            thread_data_array[i].k = k;
            pthread_create(&threads[i], NULL, parallel_portion, (void*)&thread_data_array[i]);
        }

        for (int i = 0; i < numThreads; i++) {
            pthread_join(threads[i], NULL);
        }
    }
    clock_gettime(CLOCK_REALTIME, &end); // Record end time
    double time_taken = end.tv_sec - start.tv_sec + (end.tv_nsec - start.tv_nsec) / 1e9; // Calculate elapsed time in seconds

    for (int i = 0; i < n; ++i) {
        permutation_matrix[i][p[i]] = 1.0;
    }
   // Printing results
    printf("Original matrix:\n");
    print(a_duplicate, n);
    printf("L matrix:\n");
    print(l, n);
    printf("U matrix:\n");
    print(u, n);
    printf("P array:\n");
    for (int i = 0; i < n; ++i) {
        printf("%d ", p[i]);
    }

//    // the multiplication function was used to verify the correctness of the LU decomposition
       matrix_multiply(permutation_matrix, a_duplicate, PA, n);
       matrix_multiply(l, u, LU, n);

    printf("\n\nPA matrix:\n");
    print(PA, n);
    printf("LU matrix:\n");
    print(LU, n);

    // Perform forward and backward substitution
    forwardSubstitution(l, b, y, n);
    backwardSubstitution(u, y, x, n);

    // Print the solution vector x
    printf("\n--- Solution Vector x ---\n");
     for (int i = 0; i < n; i++) {
         printf("x[%d] = %f\n", i, x[i]);
     }

    fprintf(stdout, "\nTime taken: %.9f seconds\n", time_taken);

    // Call the freeUpMemory function to deallocate memory
    freeUpMemory(&a, &a_duplicate, &u, &l, &permutation_matrix, &PA, &LU, &p);
    return 0;
}
