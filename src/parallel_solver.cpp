// Compile with: g++ -lpthread src2/parallel_solver.cpp -o bin/parallel_solver -std=c++11

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <pthread.h>
#include <vector>

using namespace std;

double **a, **a_duplicate, **u, **l, **permutation_matrix, **PA, **LU;
double *b, *x, *y;  // b is the last line of the matrix file, x and y are solution vectors
int n;
int numThreads;     // Number of threads to be used

timespec start, end_time; // Timing variables

struct thread_data {
    int id;         // thread ID
    int k;          // current step of the decomposition
};

void forwardSubstitution(double** L, double* b, double* y, int n) {
    for (int i = 0; i < n; i++) {
        y[i] = b[i];
        for (int k = 0; k < i; k++)
            y[i] -= L[i][k] * y[k];
        y[i] = y[i] / L[i][i];
    }
}

void backwardSubstitution(double** U, double* y, double* x, int n) {
    for (int i = n - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < n; j++) {
            x[i] -= U[i][j] * x[j];
        }
        x[i] /= U[i][i];
    }
}

void print(double** matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << matrix[i][j] << "  ";
        }
        cout << endl;
    }
    cout << endl;
}

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

void* parallel_portion(void* thread_data) {
    struct thread_data* my_data;
    my_data = (struct thread_data*) thread_data;
    int id = my_data->id;
    int k = my_data->k;

    int interation_per_thread = n - 1 - k;
    int start = (k + 1) + id * interation_per_thread / numThreads;
    int end = (k + 1) + (id + 1) * interation_per_thread / numThreads < n ? (k + 1) + (id + 1) * interation_per_thread / numThreads : n;
    for (int i = start; i < end; i++) {
        for (int j = k + 1; j < n; j++) {
            a[i][j] -= l[i][k] * u[k][j];
        }
    }
    pthread_exit(NULL);
}

void readMatrixFromFile(const char* filename) {
    FILE* file = fopen(filename, "r");
    if (file == NULL) {
        perror("Error opening file");
        exit(EXIT_FAILURE);
    }

    fscanf(file, "%d", &n); // Read the size of the matrix from the first line

    a = new double*[n];
    a_duplicate = new double*[n];
    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        a_duplicate[i] = new double[n];
        for (int j = 0; j < n; j++) {
            fscanf(file, "%lf", &a[i][j]);
            a_duplicate[i][j] = a[i][j]; // Copy the value to a_duplicate
        }
    }

    // Assuming the vector b is in the last line after the matrix
    b = new double[n];
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
        cerr << "This is a singular matrix, LU Decomposition is not possible" << endl;
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
        delete[] (*a)[i];
        delete[] (*a_duplicate)[i];
        delete[] (*u)[i];
        delete[] (*l)[i];
        delete[] (*permutation_matrix)[i];
        delete[] (*PA)[i];
        delete[] (*LU)[i];
    }
    delete[] *a;
    delete[] *a_duplicate;
    delete[] *u;
    delete[] *l;
    delete[] *permutation_matrix;
    delete[] *PA;
    delete[] *LU;
    delete[] *p;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        cout << "Usage: " << argv[0] << " <file path> <number of numThreads>" << endl;
        return 1;
    }
    clock_gettime(CLOCK_REALTIME, &start); // Record start time

    readMatrixFromFile(argv[1]);
    numThreads = atoi(argv[2]);

    u = new double*[n];
    l = new double*[n];
    permutation_matrix = new double*[n];
    PA = new double*[n];
    LU = new double*[n];
    x = new double[n];
    y = new double[n];

    for (int i = 0; i < n; i++) {
        u[i] = new double[n];
        l[i] = new double[n];
        permutation_matrix[i] = new double[n]();
        PA[i] = new double[n];
        LU[i] = new double[n];
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

    int* p = new int[n];
    for (int i = 0; i < n; i++) {
        p[i] = i;
    }

    pthread_t* threads = new pthread_t[numThreads];
    thread_data* thread_data_array = new thread_data[numThreads];

    for (int k = 0; k < n; k++) {
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
    clock_gettime(CLOCK_REALTIME, &end_time); // Record end time
    double time_taken = end_time.tv_sec - start.tv_sec + (end_time.tv_nsec - start.tv_nsec) / 1e9; // Calculate elapsed time in seconds

    for (int i = 0; i < n; ++i) {
        permutation_matrix[i][p[i]] = 1.0;
    }

    cout << "Original matrix:" << endl;
    print(a_duplicate, n);
    cout << "L matrix:" << endl;
    print(l, n);
    cout << "U matrix:" << endl;
    print(u, n);
    cout << "P array:" << endl;
    for (int i = 0; i < n; ++i) {
        cout << p[i] << " ";
    }
    cout << endl;

    matrix_multiply(permutation_matrix, a_duplicate, PA, n);
    matrix_multiply(l, u, LU, n);

    cout << "\n\nPA matrix:" << endl;
    print(PA, n);
    cout << "LU matrix:" << endl;
    print(LU, n);
   

    // Additional logic for applying permutation to b, performing forward and backward substitution
    double* b_prime = new double[n];
    for (int i = 0; i < n; i++) {
        b_prime[i] = b[p[i]];
    }

    forwardSubstitution(l, b_prime, y, n);
    backwardSubstitution(u, y, x, n);

    cout << "Solution vector x:" << endl;
    for (int i = 0; i < n; i++) {
        cout << "x[" << i << "] = " << x[i] << endl;
    }

    cout << "\nTime taken: " << time_taken << " seconds" << endl;

    freeUpMemory(&a, &a_duplicate, &u, &l, &permutation_matrix, &PA, &LU, &p);
    delete[] threads;
    delete[] thread_data_array;

    return 0;
}
