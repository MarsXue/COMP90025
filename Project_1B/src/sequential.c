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
#include <strings.h>

long int max(long int x, long int y) {
   return (x > y) ? x : y;
}

/* (No longer from the URL given in line 2) */
long int knapSack(long int C, long int w[], long int v[], int n) {
    int i;
    unsigned char used [n], solution [n];
    int done = 0;
    long int wt;

    bzero (used, sizeof (used));
    long int max_value = 0;
    long int weight_of_max = 0;
    long int weight = 0;
    long int value = 0;
    while (!done) {
        int carry = 1;
        done = 1;
        for (i = 0; i < n; i++)
            if (!used[i]) {
                used[i] = 1;
                weight += w[i];
                value += v[i];
                done= 0;
                break;
            } else {
                used[i] = 0;
                weight -= w[i];
                value -= v[i];
            }
        if (weight <= C && value > max_value) {
            max_value = value;
            weight_of_max = weight;
            bcopy (used, solution, sizeof (used));
        }
    }
    return max_value;
}
