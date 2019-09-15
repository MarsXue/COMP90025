/* Knapsack calculation based on that of */
/* https://www.tutorialspoint.com/cplusplus-program-to-solve-knapsack-problem-using-dynamic-programming */

#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>
#include <mpi.h>

long int knapSack(long int C, long int w[], long int v[], int n);

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(int argc, char *argv[]) {
    long int C;    /* capacity of backpack */
    int n;    /* number of items */
    int i;    /* loop counter */

    MPI_Init (&argc, &argv);
    int rank;
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        scanf ("%ld", &C);
        scanf ("%d", &n);
    }

    long int v[n], w[n];        /* value, weight */

    if (rank == 0) {
        for (i = 0; i < n; i++) {
            scanf ("%ld %ld", &v[i], &w[i]);
        }

    }

    uint64_t start = GetTimeStamp ();
    long int ks = knapSack(C, w, v, n); 

    if (rank == 0) {
        printf ("knapsack occupancy %ld\n", ks);
        printf ("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));
    }

    MPI_Finalize ();

    return 0;
}

/* PLACE YOUR CHANGES BELOW HERE */

/******************************************************************************/

/*
 * Created for COMP90025 Parallel and Multicore Computing - Project 1B, 2019
 * by Hanbin Li <hanbinl1>, Wenqing Xue <wenqingx>
 * 
 * The project is about solving the 0-1 knapsack problem in a parallel manner,
 * which is a classic NP-hard resource allocation problem. Our implementations
 * use two differetn approaches: dynamic programming with OpenMP and parallel
 * dynamic programming. By contrast, dynamic programming with OpenMP performs
 * better than the other. The detailed analysis of two methodologies will be
 * discussed in our report.
 */

#include <omp.h>
#include <stdlib.h>
#include <strings.h>

long int max(long int x, long int y) {
    return (x > y) ? x : y;
}

// MPI parallel dynamic programming approach
// For each node, they will compute [rank % size] th rows
long int parallel_dp(long int C, long int w[], long int v[], int n) {
    // broadcast the arrays of values and weights first
    MPI_Bcast(v, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    // variables for MPI send and receive opertaion
    MPI_Status status;
    MPI_Request request;

    int i;
    long int j, k;
    long int cost = 0;
    long int max_value = 0;

    // 2D array memeory allocation
    long int *K[n]; 
    for (i=0; i<n; i++) {
        K[i] = (long int *)malloc((C+1) * sizeof(long int)); 
    }

    // MPI rank, size
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // for each item from 0 to n-1
    for (i=0; i<n; i++) {
        // for each possible capacity from 0 to C
        // each rank computes its own rows
        for (j=rank; j<=C; j+=size) {
            // find the best cost value for this item
            if (j < w[i]) {
                cost = 0;
            } else if (i == 0 && j >= w[i]) {
                cost = v[i];
            } else {
                MPI_Recv(&cost, 1, MPI_LONG, (j-w[i])%size, i-1, MPI_COMM_WORLD, &status);
                cost += v[i];
            }
            // assign the maximum cost value in current situation
            // which is the same as in sequential dynamic programming
            if (i == 0) {
                K[i][j] = max(cost, 0);
            } else {
                K[i][j] = max(cost, K[i-1][j]);
            }
            // MPI_Isend is an asynchronization operation
            // send to all nodes that may need this updated value
            if (w[i+1]+j <= C) {
                MPI_Isend(&K[i][j], 1, MPI_INT, (w[i+1]+j) % size, i, MPI_COMM_WORLD, &request);
            }
        }
    }

    // make sure all nodes are completed at this point
    MPI_Barrier(MPI_COMM_WORLD);

    // only the node computes the last row will handle the final result
    // which is the right bottom corner of the 2D array
    if (rank == C % size) {
        max_value = K[n-1][C];
    }

    // free the memory allocation
    for (i=0; i<n; i++) {
        free(K[i]); 
    }

    return max_value;
}

/* (No longer from the URL given in line 2) */
long int knapSack(long int C, long int w[], long int v[], int n) {
    // MPI rank
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int i;
    long int j;

    // 2D array memeory allocation
    long int *K[n+1]; 
    for (i=0; i<=n; i++) {
        K[i] = (long int *)malloc((C+1) * sizeof(long int)); 
    }

    // iterate each item
    for (i=0; i<=n; i++) {
        // no dependency in this for loop
        // so use OpenMP for optimisation
        #pragma omp parallel for
        for (j = 0; j <= C; j++) {
            if (i == 0 || j == 0) {
                K[i][j] = 0;
            } else if (w[i-1] <= j) {
                K[i][j] = max(v[i-1] + K[i-1][j-w[i-1]], K[i-1][j]);
            } else {
                K[i][j] = K[i-1][j];
            }
        }
    }

    // make sure all nodes are completed at this point
    MPI_Barrier(MPI_COMM_WORLD);

    int max_value = 0;
    // only root node extracts the final result
    if (rank == 0) {
        max_value = K[n][C];
    }

    // free the memory allocation
    for (i=0; i<=n; i++) {
        free(K[i]); 
    }

    return max_value;
}

/* mpicc -O3 -std=c99 -fopenmp wenqingx-knapsack.c -o wenqingx-knapsack */