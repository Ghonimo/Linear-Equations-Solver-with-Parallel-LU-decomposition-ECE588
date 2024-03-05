
//----------------------------------------------------------------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <pthread.h>
#include <sys/time.h>

// Struct to replace std::pair<int, int> for thread function arguments
typedef struct {
    int id;
    int k;
} ThreadArg;

double **a, **a_copy, **u, **l, **p_matrix, **PA, **LU;
int n;
int cores;

// Function prototypes
void *assignParallel(void *arg);
void printMatrix(double **matrix, int n);
void multiplyMatrices(double **m1, double **m2, double **res, int n);
void freeMatrix(double **matrix, int n);

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Usage: %s <file_path>\n", argv[0]);
        return EXIT_FAILURE;
    }

    const int cores = 4; // Hardcoded number of cores

    // Open the file
    FILE *file = fopen(argv[1], "r");
    if (!file) {
        printf("Error opening file.\n");
        return EXIT_FAILURE;
    }

    // Read matrix size from the first line of the file
    fscanf(file, "%d", &n);

    // Allocate memory for matrices
    a = (double **)malloc(n * sizeof(double *));
    a_copy = (double **)malloc(n * sizeof(double *));
    u = (double **)malloc(n * sizeof(double *));
    l = (double **)malloc(n * sizeof(double *));
    p_matrix = (double **)malloc(n * sizeof(double *));
    PA = (double **)malloc(n * sizeof(double *));
    LU = (double **)malloc(n * sizeof(double *));

    for (int i = 0; i < n; i++) {
        a[i] = (double *)malloc(n * sizeof(double));
        a_copy[i] = (double *)malloc(n * sizeof(double));
        u[i] = (double *)malloc(n * sizeof(double));
        l[i] = (double *)malloc(n * sizeof(double));
        p_matrix[i] = (double *)malloc(n * sizeof(double));
        PA[i] = (double *)malloc(n * sizeof(double));
        LU[i] = (double *)malloc(n * sizeof(double));
    }

    // Initialize matrices and vectors
    srand48(time(NULL));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fscanf(file, "%lf", &a[i][j]) != 1) {
                printf("Error reading matrix from file.\n");
                fclose(file);
                return EXIT_FAILURE;
            }
            a_copy[i][j] = a[i][j];
            u[i][j] = (j >= i) ? drand48() : 0.0;
            l[i][j] = (j < i) ? drand48() : (j == i) ? 1.0 : 0.0;
            p_matrix[i][j] = 0.0;
        }
    }
    printMatrix(a,n);
    fclose(file);
    // Start timing
    struct timeval start, end;
    gettimeofday(&start, NULL);

///////// ---------------------- Parallel LU Decomposition ----------------------
    for (int k = 0; k < n; k++) {
        // Find the pivot: the maximum element in the k-th column starting from row k
        int kdash = k; // Index of the maximal element
        double maxi = fabs(a[k][k]);
        for (int i = k + 1; i < n; i++) {
            if (fabs(a[i][k]) > maxi) {
                maxi = fabs(a[i][k]);
                kdash = i;
            }
        }

        // Check for a singular matrix
        if (maxi == 0.0) {
            printf("Singular matrix.\n");
            exit(EXIT_FAILURE);
        }

        // Swap rows in a and a_copy
        double *temp = a[k];
        a[k] = a[kdash];
        a[kdash] = temp;

        temp = a_copy[k];
        a_copy[k] = a_copy[kdash];
        a_copy[kdash] = temp;

        // Swap rows in l for elements below diagonal (k)
        for (int i = 0; i < k; i++) {
            double tempL = l[k][i];
            l[k][i] = l[kdash][i];
            l[kdash][i] = tempL;
        }

        // Perform the division step for the l matrix
        u[k][k] = a[k][k];
        for (int i = k + 1; i < n; i++) {
            l[i][k] = a[i][k] / u[k][k];
            u[k][i] = a[k][i];
        }

        // Update the remaining elements in parallel
        pthread_t threads[cores];
        ThreadArg args[cores];
        for (int i = 0; i < cores; i++) {
            args[i].id = i;
            args[i].k = k;
            pthread_create(&threads[i], NULL, assignParallel, (void *)&args[i]);
        }

        for (int i = 0; i < cores; i++) {
            pthread_join(threads[i], NULL);
        }
    }
///////// ---------------------- Parallel LU Decomposition ----------------------
        printf("**********\n");
        printMatrix(l,n);
        printf("**********\n");
        printMatrix(u,n);
        printf("**********\n");
        multiplyMatrices(l, u, LU, n);
        printf("**********\n");
        printMatrix(LU,n);
    // End timing
    gettimeofday(&end, NULL);
    double time_spent = (end.tv_sec - start.tv_sec) + (end.tv_usec - start.tv_usec) / 1000000.0;
    printf("PThread Time: %f\n", time_spent);

    // Cleanup
    freeMatrix(a, n);
    freeMatrix(a_copy, n);
    freeMatrix(u, n);
    freeMatrix(l, n);
    freeMatrix(p_matrix, n);
    freeMatrix(PA, n);
    freeMatrix(LU, n);

    return 0;
}

void *assignParallel(void *arg) {
    ThreadArg *args = (ThreadArg *)arg;
    int id = args->id;
    int k = args->k;
    int total = n - 1 - k;
    int start = (k + 1) + id * total / cores;
    int end = (k + 1) + (id + 1) * total / cores;
    for (int i = start; i < end; i++) {
        for (int j = k + 1; j < n; j++) {
            a[i][j] -= l[i][k] * u[k][j];
        }
    }
    return NULL;
}

void printMatrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void multiplyMatrices(double **m1, double **m2, double **res, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            res[i][j] = 0;
            for (int k = 0; k < n; k++) {
                res[i][j] += m1[i][k] * m2[k][j];
            }
        }
    }
}

void freeMatrix(double **matrix, int n) {
    for (int i = 0; i < n; i++) {
        free(matrix[i]);
    }
    free(matrix);
}
