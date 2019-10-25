/*
 * Created for COMP90025 Parallel and Multicore Computing - Project 2, 2019
 * by Hanbin Li <hanbinl1>, Wenqing Xue <wenqingx>
 * 
 * The project is to simulate the N-Body problem in a parallel manner,
 * which is a well-known topic in physics and astronomy area.
 * 
 * This algorithm is implemented using openMP + MPI with a complexity of O(N^2).
 */

#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

/******************************************************************************/

typedef double vector[2];

int rank, size;
vector *v = NULL;

MPI_Datatype MPI_VECTOR;

// function is based on previous project
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

// read the input from the command line arguments
void read_input(char *argv[], int *N, int *T, double *G, double *delta_t) {
    // root node reads the input
    if (rank == 0) {
        *N = strtol(argv[1], NULL, 10);     // the number of particles
        *T = strtol(argv[2], NULL, 10);     // the number of time steps
        *G = strtod(argv[3], NULL);         // the constant gravity force
        *delta_t = strtod(argv[4], NULL);   // the value of delta time
    }

    // broadcasts the data all other processes of the communicator
    MPI_Bcast(N, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(T, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(G, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(delta_t, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
}

// read the input data file from given file name
void read_file(char *file_name, double m[], vector p[], vector recv_v[], int N, int chunk) {
    if (rank == 0) {
        // open the file
        FILE *file = fopen(file_name, "r");

        // check the input file existed
        if (!file) {
            printf("No input file <%s>\n", file_name);
            MPI_Finalize();
            exit(1);
        }

        // scan the mass, position x, position y, velocity x, velocity y
        for (int i=0; i<N; i++) {
            fscanf(file, "%lf %lf %lf %lf %lf", &m[i],
                                                &p[i][0], &p[i][1],
                                                &v[i][0], &v[i][1]);
        }
        // close the file
        fclose(file);
    }
    // broadcasts the data all other processes of the communicator
    MPI_Bcast(m, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(p, N, MPI_VECTOR, 0, MPI_COMM_WORLD);
    MPI_Scatter(v, chunk, MPI_VECTOR, recv_v, chunk, MPI_VECTOR, 0, MPI_COMM_WORLD);
}


// compute the acceleration of the particles
void compute_force(int idx, double m[], vector p[], vector f[], int N, double G, int chunk) {

    vector tmp;
    double distance, factor;

    // the index of current particle
    int curr = rank * chunk + idx;

    // initialize to zero for summation
    f[idx][0] = 0.0;
    f[idx][1] = 0.0;

    for (int i = 0; i < N; i++) {
        if (i != curr) {
            // difference of x coordinate
            tmp[0] = p[curr][0] - p[i][0];
            // difference of y coordinate
            tmp[1] = p[curr][1] - p[i][1];
            // square root of the sum of differences' power
            distance = sqrt(pow(tmp[0], 2) + pow(tmp[1], 2));
            // partial equation in the report
            // G * m_i * m_j / (p_j - p_i)^3
            factor = G * m[curr] * m[i] / pow(distance, 3);
            // G * m_i * m_j * (p_j - p_i) / (p_j - p_i)^3
            tmp[0] *= factor;
            tmp[1] *= factor;
            // sum up
            f[idx][0] += tmp[0];
            f[idx][1] += tmp[1];
        }
    }
}

// update the positions and velocities of particles
void update_particles(int idx, double m[], vector f[], vector recv_p[], vector recv_v[], int N, int chunk, double delta_t) {
    int curr = rank * chunk + idx;
    double factor = delta_t / m[curr];

    // update position x
    recv_p[idx][0] += delta_t * recv_v[idx][0];
    // update position y
    recv_p[idx][1] += delta_t * recv_v[idx][1];
    // update velocity x
    recv_v[idx][0] += factor * f[idx][0];
    // update velocity y
    recv_v[idx][1] += factor * f[idx][1];
}

// print the output data of simulated particles
void print_output(double m[], vector p[], vector recv_v[], int N, int chunk) {
    MPI_Gather(recv_v, chunk, MPI_VECTOR, v, chunk, MPI_VECTOR, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        for (int i = 0; i < N; i++) {
            printf("%lf %lf %lf %lf %lf\n", m[i],
                                            p[i][0], p[i][1],
                                            v[i][0], v[i][1]);
        }
    }
}

// main program
int main(int argc, char *argv[]) {
    int N;                      // the number of particles
    int T;                      // the number of time steps
    double G;                   // the input gravity force
    double delta_t;             // the number of delta time

    int chunk;                  // the value of chunk size
    char *file_name;            // the input data file name

    double *m;                  // the masses of particles
    vector *p;                  // the positions of paticles
    vector *f;                  // the forces of particles
    vector *recv_p;             // the received position
    vector *recv_v;             // the received velocity

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

    // chunk is the number of particles divides the number of processor
    chunk = N / size;

    file_name = argv[5];

    // malloc the memory
    m = malloc(N * sizeof(double));
    p = malloc(N * sizeof(vector));
    f = malloc(chunk * sizeof(vector));
    recv_p = p + rank * chunk;
    recv_v = malloc(chunk * sizeof(vector));

    if (rank == 0) {
        v = malloc(N * sizeof(vector));
    }

    MPI_Type_contiguous(2, MPI_DOUBLE, &MPI_VECTOR);
    MPI_Type_commit(&MPI_VECTOR);

    // read the input file
    read_file(file_name, m, p, recv_v, N, chunk);

    // at each time step
    for (int step = 0; step < T; step++) {
        // compute the force
        #pragma omp for schedule(static, N / size)
        for (int i = 0; i < chunk; i++) {
            compute_force(i, m, p, f, N, G, chunk);
        }
        // update the particles
        #pragma omp for schedule(static, N / size)
        for (int i = 0; i < chunk; i++) {
            update_particles(i, m, f, recv_p, recv_v, N, chunk, delta_t);
        }
        // gathers data from all tasks and distribute the combined data to all tasks
        MPI_Allgather(MPI_IN_PLACE, chunk, MPI_VECTOR, p, chunk, MPI_VECTOR, MPI_COMM_WORLD);
    }

    // print the output simulated data
    print_output(m, p, recv_v, N, chunk);

    if (rank == 0) {
        printf ("Time: %llu us\n", (uint64_t) (GetTimeStamp() - start));
        free(v);
    }

    // free the memory
    free(m);
    free(p);
    free(f);
    free(recv_v);

    // free the MPI type
    MPI_Type_free(&MPI_VECTOR);
    MPI_Finalize();

    return 0;
}

/* mpicc -O3 -std=c99 -fopenmp N2_code.c -o N2_code -lm */
