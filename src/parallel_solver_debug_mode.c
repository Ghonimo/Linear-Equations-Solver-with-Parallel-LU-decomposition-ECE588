// cc -lpthread src/parallel_solver.c -o bin/parallel_solver -std=gnu11

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>

// Function prototypes --------------------------
void forwardSubstitution(double** L, double* b, double* y, int n);
void backwardSubstitution(double** U, double* y, double* x, int n);
void print(double** matrix, int n);
void matrix_multiply(double** matrix1, double** matrix2, double** product, int n);
void readMatrixFromFile(const char* filename);
void performPartialPivoting(double** a, double** l, int* p, int k);
void freeUpMemory(double*** a, double*** a_duplicate, double*** u, double*** l, double*** permutation_matrix, double*** PA, double*** LU, int** p);
void *parallel_portion(void* thread_data);
double ludecomp_verify(double** P, double** A, double** L, double** U, double** residual, int n);
// ----------------------------------------------

double **a, **a_duplicate, **u, **l, **permutation_matrix, **PA, **LU;
double *b, *x, *y;  // b is the last line of the matrix file, x and y are solution vectors
int n;
int numThreads;     // Number of threads to be used

struct timespec start, end; // Timing variables
struct thread_data {
    int id;         // thread ID
    int k;          // current step of the decomposition
};

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

    // Allocate memory for the residual matrix 'residual'
    double** residual = (double**)malloc(n * sizeof(double*));
    for (int i = 0; i < n; i++) {
        residual[i] = (double*)malloc(n * sizeof(double));
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
    printf("Completed LU Decomposition with %d threads for %dx%d matrix.\n", numThreads, n, n);
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

   // the multiplication function was used to verify the correctness of the LU decomposition
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

    double norm = ludecomp_verify(permutation_matrix, a_duplicate, l, u, residual, n);
    printf("\n Frobenius norm of the residual matrix: %f\n", norm);

    fprintf(stdout, "\nTime taken: %.9f seconds\n", time_taken);

    // Free up memory for the residual matrix 'residual'
    for (int i = 0; i < n; i++) {
        free(residual[i]);
    }
    free(residual);

    // Call the freeUpMemory function to deallocate memory
    freeUpMemory(&a, &a_duplicate, &u, &l, &permutation_matrix, &PA, &LU, &p);
    return 0;
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

    // Read the size of the matrix from the first line
    fscanf(file, "%d", &n); 

    // Reading the matrix from the file
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
    
    // read vector b from the last line after the matrix
    b = (double*)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++) {
        fscanf(file, "%lf", &b[i]);
    }

    fclose(file);
}

/**
 * Executes LU decomposition steps in parallel by updating matrix A.
 *
 * This function is designed to run as a separate thread, where each thread
 * works on a different portion of matrix A during the LU decomposition process.
 * It updates the elements of A based on the current decomposition step, represented
 * by the index k. The work is divided among threads so that each thread updates
 * a unique set of rows, ensuring parallel computation without data races on distinct rows.
 *
 * Parameters:
 *  - thread_data: Pointer to struct containing thread ID and the current decomposition step (k).
 *                 The struct should have fields 'id' for the thread ID and 'k' for the step index.
 *
 * Process:
 *  - Retrieves its assigned portion of the matrix based on 'k' and its thread ID ('id').
 *  - Iterates over its portion, updating elements in matrix A using values from L and U matrices.
 *    The formula applied is A[i][j] -= L[i][k] * U[k][j] for each element in the assigned range.
 *  - Prints updates to A[i][j] for transparency and debugging purposes.
 *
 * This parallel portion aims to expedite the LU decomposition process by leveraging
 * multi-threading, where each thread contributes to computing parts of the LU factorization.
 * 
 * Note:
 *  - 'n' is the global variable representing the size of the square matrix.
 *  - 'numThreads' is the global variable indicating how many threads are participating
 *    in the decomposition process.
 */
void* parallel_portion(void* thread_data) {
    struct thread_data* my_data;
    my_data = (struct thread_data*) thread_data;
    int id = my_data->id;
    int k = my_data->k;

    printf("Thread %d: Computing step k: %d\n", id, k);
    int rows_per_step = n - 1 - k;
    int start = (k + 1) + id * rows_per_step / numThreads;
    int end = (k + 1) + (id + 1) * rows_per_step / numThreads < n ? (k + 1) + (id + 1) * rows_per_step / numThreads : n;
    printf("Thread %d: row from: %d to: %d\n", id, start, end);
    for (int i = start; i < end; i++) {
        for (int j = k + 1; j < n; j++) {
            a[i][j] -= l[i][k] * u[k][j]; // Update element A[i][j] for LU decomposition.
            printf("Thread %d: A[%d][%d] = %f\n", id, i, j, a[i][j]);
        }
    }
    pthread_exit(NULL); // Terminate the thread after completing its task.
}

/**
 * Performs forward substitution to solve the equation Ly = b for y.
 * where L is a lower triangular matrix. This process takes advantage of the
 * triangular structure to solve for each component of y sequentially, starting
 * from the top row. The algorithm subtracts the known contributions of previously
 * calculated y values from the current b value, then divides by the diagonal coefficient
 * of L to solve for the current y. This method ensures a solution with a linear time
 * complexity relative to the number of equations, which is highly efficient for
 * lower triangular matrices.
 *
 * Parameters:
 *  - L: A double pointer to the lower triangular matrix L, where L[i][j] = 0 for all j > i.
 *       The matrix is assumed to be non-singular (no zero elements on the diagonal).
 *  - b: A pointer to the vector b, representing the right-hand side of the equation Ly = b.
 *  - y: A pointer to the solution vector y, where the solution will be stored.
 *       The memory for y should be allocated before calling this function.
 *  - n: The size of the matrix L and vectors b and y, representing the number of linear equations.
 *
 * The function directly modifies the array pointed to by y to contain the solution of Ly = b.
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
 * Performs backward substitution to solve the equation Ux = y for x.
 * where U is an upper triangular matrix. This algorithm utilizes the upper triangular
 * structure by starting the solution process from the bottom row and moving upwards,
 * solving for each component of x in a sequential manner. It subtracts the contributions
 * of previously solved x values from the current y value and then divides by the diagonal
 * coefficient of U to solve for the current x. This backward approach ensures that the
 * solution can be found efficiently, leveraging the matrix's structure to reduce computation.
 *
 * Parameters:
 *  - U: A double pointer to the upper triangular matrix U, where U[i][j] = 0 for all j < i.
 *       The matrix is assumed to be non-singular (no zero elements on the diagonal).
 *  - y: A pointer to the vector y, representing the right-hand side of the equation Ux = y.
 *  - x: A pointer to the solution vector x, where the solution will be stored.
 *       The memory for x should be allocated before calling this function.
 *  - n: The size of the matrix U and vectors y and x, representing the number of linear equations.
 *
 * The function modifies the array pointed to by x to contain the solution of Ux = y.
 * The solution is computed in reverse order, starting from the last equation and
 * moving up to the first, effectively utilizing the upper triangular form of U
 * to simplify the calculation process.
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
 * Prints a square matrix.
 * 
 * This function iterates through a square matrix, printing each element
 * formatted as a floating point number. Rows are printed on new lines.
 *
 * @param matrix A double pointer to the square matrix to be printed.
 * @param n The dimension of the square matrix.
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
 * Multiplies two square matrices.
 * 
 * This function computes the product of two square matrices (matrix1 and matrix2)
 * and stores the result in 'product'. The multiplication follows the standard rule
 * of matrix multiplication, iterating through rows of the first matrix and columns
 * of the second matrix.
 *
 * @param matrix1 A double pointer to the first square matrix.
 * @param matrix2 A double pointer to the second square matrix.
 * @param product A double pointer to the square matrix where the result is stored.
 * @param n The dimension of the square matrices.
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
 * Performs partial pivoting on matrices for LU decomposition.
 *
 * This function searches for the maximum absolute value in column k below or at
 * the diagonal and swaps the current row with the row containing this maximum
 * to maintain numerical stability in LU decomposition. It also updates the
 * permutation vector to reflect the row swaps.
 *
 * @param a A double pointer to the matrix being decomposed.
 * @param l A double pointer to the lower triangular matrix being formed.
 * @param p An integer pointer to the permutation vector that records the row swaps.
 * @param k The current column (and row) index at which the pivoting is performed.
 */
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

/**
 * Frees dynamically allocated memory for matrices and vectors.
 *
 * This function releases the memory allocated for matrices and vectors used
 * during LU decomposition and solving linear systems, preventing memory leaks.
 *
 * @param a A pointer to a double pointer to the original matrix.
 * @param a_duplicate A pointer to a double pointer to the duplicate of the original matrix.
 * @param u A pointer to a double pointer to the upper triangular matrix.
 * @param l A pointer to a double pointer to the lower triangular matrix.
 * @param permutation_matrix A pointer to a double pointer to the permutation matrix.
 * @param PA A pointer to a double pointer to the product of permutation and original matrices.
 * @param LU A pointer to a double pointer to the product of L and U matrices.
 * @param p A pointer to an integer pointer to the permutation vector.
 */
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

/**
 * Calculates the residual norm to verify the correctness of LU Decomposition.
 *
 * This function computes the residual matrix by subtracting the product of
 * L and U matrices from the product of the permutation matrix P and the original
 * matrix A. The residual matrix reflects the difference between the left-hand
 * and right-hand sides of the PA=LU equation. The function then calculates and
 * returns the Frobenius norm of the residual matrix as a measure of the decomposition's
 * accuracy. A smaller norm suggests a more accurate decomposition.
 *
 * Parameters:
 *  - P: A double pointer to the permutation matrix.
 *  - A: A double pointer to the original matrix.
 *  - L: A double pointer to the lower triangular matrix.
 *  - U: A double pointer to the upper triangular matrix.
 *  - residual: A double pointer to the residual matrix where the residual will be stored.
 *  - n: The size of the matrices, representing an n x n system.
 *
 * Returns:
 *  - The Frobenius norm of the residual matrix, which quantifies the accuracy of the LU Decomposition.
 */
double ludecomp_verify(double** P, double** A, double** L, double** U, double** residual, int n) {
    // Compute the residual matrix: (PA - LU)
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            residual[i][j] = 0;
            for (int k = 0; k < n; k++) {
                residual[i][j] += P[i][k] * A[k][j];          // Adding product of P and A matrices
                residual[i][j] -= L[i][k] * U[k][j];          // Subtracting product of L and U matrices
            }
        }
    }

    // Calculate the Frobenius norm of the residual matrix
    double norm = 0;
    for (int i = 0; i < n; i++) {
        double column_norm = 0;
        for (int j = 0; j < n; j++) {
            column_norm += residual[j][i] * residual[j][i]; // Summing squares of each element in a column
        }
        norm += sqrt(column_norm); // Adding the square root of the summed squares of each column
    }
    return norm;                    // Return the Frobenius norm as the measure of accuracy
}
