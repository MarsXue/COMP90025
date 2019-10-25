/*
 * Created for COMP90025 Parallel and Multicore Computing - Project 2, 2019
 * by Hanbin Li <hanbinl1>, Wenqing Xue <wenqingx>
 * 
 * The project is to simulate the N-Body problem in a parallel manner,
 * which is a well-known topic in physics and astronomy area.
 *
 * This algorithm is implemented using Barnes Hut algorithm with OpenMP+MPI.
 * This approach has a complexity of O(N*log(N)).
 * There is an additional file "tree.c" used for BH algorithm.
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include "tree.c"

/******************************************************************************/

const double treeratio = 1;

int rank, size;

double max(double a, double b) {
    return (a > b ? a : b);
}

double min(double a, double b) {
    return (a < b ? a : b);
}

// function is based on previous project
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

// read the input from the command line arguments
void read_input(char *argv[], int *N, int *T, double *G, double *delta_t){
    if (rank == 0) {
        *N = strtol(argv[1], NULL, 10);
        *T = strtol(argv[2], NULL, 10);
        *G = strtod(argv[3], NULL);   
        *delta_t = strtod(argv[4], NULL);
    }
    // broadcasts the data all other processes of the communicator
    MPI_Bcast(N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(T, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(G, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(delta_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}


// read the initial conditions into bodies struct
body_t *read_file(char *file_name, body_t *bodies, int N){
    
    FILE * file = fopen(file_name, "r");
    if (!file) {
        printf("No input file <%s>\n", file_name);
        MPI_Finalize();
        exit(1);
    }
    
    // malloc the memory
    bodies = malloc(N * sizeof(body_t));
    
    // scan the mass, position x, position y, velocity x, velocity y
    for (int i = 0; i < N; i++) {
        fscanf(file, "%lf %lf %lf %lf %lf", &((bodies)[i].m), 
                                            &((bodies)[i].px), &((bodies)[i].py),
                                            &((bodies)[i].vx), &((bodies)[i].vy));
    }
    // close the file
    fclose(file);
    return bodies;
}

// print the output data of simulated particles
void print_output(body_t *bodies, int N) {
    if (rank == 0) {
        for (int i = 0; i < N;i++) {
            printf("%lf %lf %lf %lf %lf\n", bodies[i].m, 
                                            bodies[i].px, bodies[i].py,
                                            bodies[i].vx, bodies[i].vy);
        }
    }
}

// main program
int main(int argc, char *argv[]) {
    int N;
    int T;
    double G;
    double delta_t;

    char *file_name;
    
    body_t *bodies = NULL;

    // initialize MPI with rank and size
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // start the timer
    uint64_t start = GetTimeStamp();

    // check the input amount is satisfied
    if (argc != 6) {
        // error message
        printf("Usage: mpirun -np <#processor> ./N2_code ");
        printf("<#particle> <time step> <gravity> <delta time> <input file>\n");
        MPI_Finalize();
        exit(0);
    }
    
    // read the input from the command line arguments
    read_input(argv, &N, &T, &G, &delta_t);

    file_name = argv[5];
    
    // read the input file
    if (rank == 0) {
        bodies = read_file(file_name, bodies, N);
    }
    
    if (rank != 0) {
        if (!(bodies = malloc(N * sizeof(body_t)))) {
            printf("Impossibile allocare memoria per %d bodies", N);
            exit(1);
        }
    }
    
    int i;
    
    //broadcast masses, positions and forces of particles to other process
    for (i = 0; i < N; i++) {
        MPI_Bcast(&((bodies)[i].m), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].px), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].py), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].fx), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(&((bodies)[i].fy), 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    }

    // at each time step
    for (int step = 0; step < T; step++) {
        
        double xmin = 0.0, xmax = 0.0;
        double ymin = 0.0, ymax = 0.0;
        
        #pragma omp for schedule(static, N / size)
        for (i = 0; i < N; i++) {
            bodies[i].fx = 0.0;
            bodies[i].fy = 0.0;
            xmin = min(xmin, bodies[i].px);
            xmax = max(xmax, bodies[i].px);
            ymin = min(ymin, bodies[i].py);
            ymax = max(ymax, bodies[i].py);
        }
        
        // construct quad-tree
        node_t *rootnode = create_node(bodies, xmin, xmax, ymin, ymax);
        
        #pragma omp for schedule(static, N / size)
        for (i = 1; i < N; i++) {
            insert_body(bodies + i, rootnode);
        }
        
        // compute force
        #pragma omp for schedule(static, N / size)
        for (i = rank; i < N; i += size) {
            tree_sum(rootnode, bodies + i, G, treeratio);
        }
        
        #pragma omp for schedule(static, N / size)
        for (i = 0; i < N; i++){
            MPI_Bcast(&(bodies[i].fx), 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
            MPI_Bcast(&(bodies[i].fy), 1, MPI_DOUBLE, i % size, MPI_COMM_WORLD);
        }
        
        //update positions and velocities
        #pragma omp for schedule(static, N / size)
        for (i = 0; i < N; i++){
            bodies[i].px += delta_t * bodies[i].vx;
            bodies[i].py += delta_t * bodies[i].vy;
            
            bodies[i].vx += bodies[i].fx * (delta_t / bodies[i].m);
            bodies[i].vy += bodies[i].fy * (delta_t / bodies[i].m);
        }
        delete_tree(rootnode);
    }
    
    // generate output
    print_output(bodies, N);
    
    
    if (rank == 0) {
        printf ("Time: %llu us\n", (uint64_t) (GetTimeStamp() - start));
    }

    // free the memory
    free(bodies);

    MPI_Finalize();

    return 0;
}

/* mpicc -O3 -std=c99 -fopenmp NlogN_code.c -o NlogN_code */
