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
    
    MPI_Bcast(&C, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
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
        printf ("Time: %llu us\n", (uint64_t) (GetTimeStamp() - start));
    }

    MPI_Finalize ();

    return 0;
}

/* PLACE YOUR CHANGES BELOW HERE */

/*
 * Created for COMP90025 Parallel and Multicore Computing - Project 1B, 2019
 * by Hanbin Li <hanbinl1>, Wenqing Xue <wenqingx>
 */
// #include <omp.h>
#include <stdlib.h>
#include <strings.h>

long int max(long int x, long int y) {
    return (x > y) ? x : y;
}

/* (No longer from the URL given in line 2) */
long int knapSack(long int C, long int w[], long int v[], int n) {
    // broadcast the arrays of values and weights first
    MPI_Bcast(v, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Bcast(w, n, MPI_LONG, 0, MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);

    MPI_Status status;
    MPI_Request request;

    int i, j, k;
    // long int cost[n][C+1];
    long int cost_i = 0;
    long int max_value = 0;

    long int *cost[n]; 
    for (i=0; i<n; i++) {
        cost[i] = (long int *)malloc((C+1) * sizeof(long int)); 
    }

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    // for each item
    for(i=0; i<n; i++){
        printf("rank %d i = %d\n", rank, i);
        // each rank computes its columns
        // #pragma omp parallel for schedule(dynamic)
        for(j=rank; j<=C; j+=size){
            
            // find the best value with the new item
            if (j - w[i] < 0) {
                cost_i = 0;
            } else if (i == 0 && j-w[i] >= 0) {
                cost_i = v[i];
            } else {
                MPI_Recv(&cost_i, 1, MPI_INT, (j-w[i]) % size, i-1, MPI_COMM_WORLD, &status);
                cost_i += v[i];
            }

            // find the maximum of current cost
            if (i == 0) {
                cost[i][j] = max(cost_i, i);
            } else {
                cost[i][j] = max(cost_i, cost[i-1][j]);
            }

            // Send to all procs that could need the new value (non blancking)
            for (k=j+1; k<=C; k++) {
                if (k-j == w[i+1]) {
                    MPI_Isend(&cost[i][j], 1, MPI_INT, k % size, i, MPI_COMM_WORLD, &request);
                }
            }
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank == C % size) {
        max_value = cost[n-1][C];
    }

    // free the malloc array
    for (i=0; i<n; i++) {
        free(cost[i]); 
    }

    return max_value;
}

/*mpicc -O3 wenqingx-knapsack.c -o wenqingx-knapsack*/