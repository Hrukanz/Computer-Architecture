#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <time.h>


/**
 * poisson.c
 * Implementation of a Poisson solver with Dirichlet boundary conditions.
 *
 * BUILDING:
 * gcc -o poisson poisson.c -lpthread
 * gcc -O3 -march=native -ffast-math -o poisson poisson.c -lpthread
 * 
 *
 * USAGE:
 * ./poisson [-n size] [-x source x-position] [-y source y-position] [-z source z-position]
 *           [-a source amplitude] [-i iterations] [-t threads] [-d]
 */

// Global flag
// Set to true when operating in debug mode to enable verbose logging
static bool debug = false;

// Structure to hold shared data among threads
typedef struct {
    int n;
    int iterations;
    double delta_sq;
    double *curr;
    double *next;
    double *source;
    pthread_barrier_t barrier;
} SharedData;

// Structure to pass arguments to worker threads
typedef struct {
    int thread_id;
    int start_k;
    int end_k;
    SharedData *shared;
} WorkerArgs;

/**
 * @brief Worker function for each thread.
 */
void* worker(void* arg) {
    WorkerArgs* args = (WorkerArgs*)arg;
    int thread_id = args->thread_id;
    int n = args->shared->n;
    int start_k = args->start_k;
    int end_k = args->end_k;
    int iterations = args->shared->iterations;
    double delta_sq = args->shared->delta_sq;
    double *curr = args->shared->curr;
    double *next = args->shared->next;
    double *source = args->shared->source;
    pthread_barrier_t *barrier = &args->shared->barrier;

    int i_plus1, i_minus1, j_plus1, j_minus1;

    for (int iter = 0; iter < iterations; ++iter) {
        for (int k = start_k; k < end_k; ++k) {
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < n; ++i) {
                    // Calculate the index into the 1D array
                    int idx = (k * n + j) * n + i;

                    // Apply boundary conditions
                    if (k == 0) {
                        // Bottom boundary (Dirichlet condition)
                        next[idx] = -1.0;
                    } else if (k == n - 1) {
                        // Top boundary (Dirichlet condition)
                        next[idx] = 1.0;
                    } else {
                        // Handle Neumann boundary conditions using mirrored indices
                        // Determine i_plus1 and i_minus1 with boundary checks
                        if (i + 1 <= n - 1) {
                            i_plus1 = i + 1;
                        } else {
                            i_plus1 = i - 1;
                        }

                        if (i - 1 >= 0) {
                            i_minus1 = i - 1;
                        } else {
                            i_minus1 = i + 1;
                        }

                        // Determine j_plus1 and j_minus1 with boundary checks
                        if (j + 1 <= n - 1) {
                            j_plus1 = j + 1;
                        } else {
                            j_plus1 = j - 1;
                        }

                        if (j - 1 >= 0) {
                            j_minus1 = j - 1;
                        } else {
                            j_minus1 = j + 1;
                        }

                        // For k+1 and k-1
                        double v_kp1;
                        if (k + 1 <= n - 1) {
                            int idx_kp1 = ((k + 1) * n + j) * n + i;
                            v_kp1 = curr[idx_kp1];
                        } else {
                            // Top boundary (Dirichlet condition)
                            v_kp1 = 1.0;
                        }

                        double v_km1;
                        if (k - 1 >= 0) {
                            int idx_km1 = ((k - 1) * n + j) * n + i;
                            v_km1 = curr[idx_km1];
                        } else {
                            // Bottom boundary (Dirichlet condition)
                            v_km1 = -1.0;
                        }

                        // Calculate the indices of the neighboring points
                        int idx_ip1 = (k * n + j) * n + i_plus1;
                        int idx_im1 = (k * n + j) * n + i_minus1;
                        int idx_jp1 = (k * n + j_plus1) * n + i;
                        int idx_jm1 = (k * n + j_minus1) * n + i;

                        // Get neighboring values
                        double v_ip1 = curr[idx_ip1];
                        double v_im1 = curr[idx_im1];
                        double v_jp1 = curr[idx_jp1];
                        double v_jm1 = curr[idx_jm1];

                        // Update the potential at the current point
                        next[idx] = (1.0 / 6.0) * (v_ip1 + v_im1 + v_jp1 + v_jm1 + v_kp1 + v_km1 - delta_sq * source[idx]);
                    }
                }
            }
        }

        // Synchronise threads
        pthread_barrier_wait(barrier);

        // Swap curr and next pointers (only in one thread)
        if (thread_id == 0) {
            double *temp = args->shared->curr;
            args->shared->curr = args->shared->next;
            args->shared->next = temp;
        }

        // Synchronise again to ensure all threads see the updated pointers
        pthread_barrier_wait(barrier);

        // Update local pointers
        curr = args->shared->curr;
        next = args->shared->next;
    }

    return NULL;
}

/**
 * @brief Solve Poissons equation for a given cube with Dirichlet boundary
 * conditions on all sides.
 *
 * @param n             The edge length of the cube. n^3 number of elements.
 * @param source        Pointer to the source term cube, a.k.a. forcing function.
 * @param iterations    Number of iterations to perform.
 * @param threads       Number of threads to use for solving.
 * @param delta         Grid spacing.
 * @return double*      Solution to Poissons equation.  Caller must free.
 */
double* poisson_mixed(int n, double *source, int iterations, int threads, float delta) {
    if (debug) {
        printf("Starting solver with:\n"
               "n = %i\n"
               "iterations = %i\n"
               "threads = %i\n"
               "delta = %f\n",
               n, iterations, threads, delta);
    }

    // Allocate buffers for computation
    double *curr = (double*)calloc(n * n * n, sizeof(double));
    double *next = (double*)calloc(n * n * n, sizeof(double));

    if (curr == NULL || next == NULL) {
        fprintf(stderr, "Error: ran out of memory when trying to allocate %i sized cube\n", n);
        exit(EXIT_FAILURE);
    }

    double delta_sq = delta * delta;

    // Initialise barrier
    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, threads);

    // Shared data among threads
    SharedData shared;
    shared.n = n;
    shared.iterations = iterations;
    shared.delta_sq = delta_sq;
    shared.curr = curr;
    shared.next = next;
    shared.source = source;
    shared.barrier = barrier;

    // Create threads
    pthread_t *thread_handles = malloc(threads * sizeof(pthread_t));
    WorkerArgs *args = malloc(threads * sizeof(WorkerArgs));

    int k_per_thread = n / threads;
    int remainder = n % threads;
    int start_k = 0;

    for (int t = 0; t < threads; ++t) {
        args[t].thread_id = t;
        args[t].shared = &shared;
        args[t].start_k = start_k;
        int k_count = k_per_thread + (t < remainder ? 1 : 0);
        args[t].end_k = start_k + k_count;
        start_k += k_count;

        if (pthread_create(&thread_handles[t], NULL, worker, (void*)&args[t]) != 0) {
            fprintf(stderr, "Error creating worker thread!\n");
            exit(EXIT_FAILURE);
        }
    }

    // Wait for all threads to complete
    for (int t = 0; t < threads; ++t) {
        pthread_join(thread_handles[t], NULL);
    }

    // Destroy barrier
    pthread_barrier_destroy(&barrier);

    free(thread_handles);
    free(args);

    if (debug) {
        printf("Finished solving.\n");
    }

    // Return the correct buffer after swapping
    if (iterations % 2 == 0) {
        // Even number of iterations, curr is the result
        free(next);
        return curr;
    } else {
        // Odd number of iterations, next is the result
        free(curr);
        return next;
    }
}

int main(int argc, char **argv)
{
    // Default settings for solver
    int iterations = 10000;
    int n = 7;
    int threads = 16;
    float delta = 1;
    int x = -1;
    int y = -1;
    int z = -1;
    double amplitude = 1.0;

    int opt;

    // parse the command line arguments
    while  ((opt = getopt(argc, argv, "h:n:i:x:y:z:a:t:d")) != -1)
    {
        switch(opt)
        {
        case 'h':
            printf("Usage: poisson [-n size] [-x source x-position] [-y source y-position] [-z source z-position] [-a source amplitude] [-i iterations] [-t threads] [-d] (for debug mode)\n");
            return EXIT_SUCCESS;
        case 'n':
            n = atoi(optarg);
            break;
        case 'i':
            iterations = atoi(optarg);
            break;
        case 'x':
            x = atoi(optarg);
            break;
        case 'y':
            y = atoi(optarg);
            break;
        case 'z':
            z = atoi(optarg);
            break;
        case 'a':
            amplitude = atof(optarg);
            break;
        case 't':
            threads = atoi(optarg);
            break;
        case 'd':
            debug = true;
            break;
        default:
            fprintf(stderr, "Usage: poisson [-n size] [-x source x-position] [-y source y-position] [-z source z-position] [-a source amplitude]  [-i iterations] [-t threads] [-d] (for debug mode)\n");
            exit(EXIT_FAILURE);
        }
    }

    // Ensure we have an odd sized cube
    if (n % 2 == 0)
    {
        fprintf(stderr, "Error: n should be an odd number!\n");
        return EXIT_FAILURE;
    }

    // Create a source term with a single point in the centre
    double *source = (double*)calloc(n * n * n, sizeof(double));
    if (source == NULL)
    {
        fprintf(stderr, "Error: failed to allocated source term (n=%i)\n", n);
        return EXIT_FAILURE;
    }

    // Default x,y, z
    if (x < 0 || x > n - 1)
        x = n / 2;
    if (y < 0 || y > n - 1)
        y = n / 2;
    if (z < 0 || z > n - 1)
        z = n / 2;

    source[(z * n + y) * n + x] = amplitude;

    // Start timer
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);

    // Calculate the resulting field with mixed boundary conditions
    double *result = poisson_mixed(n, source, iterations, threads, delta);

    // End timer
    clock_gettime(CLOCK_MONOTONIC, &end);

    // Calculate elapsed time
    double elapsed_time = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec)/1e9;

    printf("Execution time: %f seconds\n", elapsed_time);

    // Print out the middle slice of the cube for validation
    for (int y = 0; y < n; ++y)
    {
        for (int x = 0; x < n; ++x)
        {
            printf("%0.5f ", result[((n / 2) * n + y) * n + x]);
        }
        printf("\n");
    }

    free(source);
    free(result);

    return EXIT_SUCCESS;
}