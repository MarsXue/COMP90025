#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define DELTA_TIME 0.01 // delta time

typedef struct {
    double m;           // mass    
    double fx, fy, fz;  // force
    double px, py, pz;  // position
    double vx, vy, vz;  // velocity
} Body;

double G;
int N, T;
Body *bodies;

// Computes the distance between two particles
double compute_distance(Body a, Body b) {
    return sqrt(pow((a.px - b.px), 2) + pow((a.py - b.py), 2) + pow((a.pz - b.pz), 2));
}

// Computes the force for each particle
void compute_force() {
    for (int i=0; i<N; i++) {
        bodies[i].fx = 0.0;
        bodies[i].fy = 0.0;
        bodies[i].fz = 0.0;
        for (int j=0; j<N; j++) {
            if (i == j) continue;
            double d = compute_distance(bodies[i], bodies[j]);
            double f = (G * bodies[i].m * bodies[j].m) / (pow(d, 2.0));
            
            bodies[i].fx += f * ((bodies[j].px - bodies[i].px) / d);
            bodies[i].fy += f * ((bodies[j].py - bodies[i].py) / d);
            bodies[i].fz += f * ((bodies[j].pz - bodies[i].pz) / d);
        }
    }
}

// Computes the velocity for each particle
void compute_velocity() {
    for (int i=0; i<N; i++) {
        bodies[i].vx += (bodies[i].fx / bodies[i].m) * DELTA_TIME;
        bodies[i].vy += (bodies[i].fy / bodies[i].m) * DELTA_TIME;
        bodies[i].vz += (bodies[i].fz / bodies[i].m) * DELTA_TIME;
    }
}

// Computes the position for each particle
void compute_position() {
    for (int i=0; i<N; i++) {
        bodies[i].px += bodies[i].vx * DELTA_TIME;
        bodies[i].py += bodies[i].vy * DELTA_TIME;
        bodies[i].pz += bodies[i].vz * DELTA_TIME;
    }
}

void resolve_collision() {
    for (int i=0; i<N-1; i++) {
        for (int j=i+1; j<N; j++) {
            if (bodies[i].px == bodies[j].px &&
                bodies[i].py == bodies[j].py &&
                bodies[i].pz == bodies[j].pz) {
                double tx = bodies[i].vx;
                double ty = bodies[i].vy;
                double tz = bodies[i].vz;
                bodies[i].vx = bodies[j].vx;
                bodies[i].vy = bodies[j].vy;
                bodies[i].vz = bodies[j].vz;
                bodies[j].vx = tx;
                bodies[j].vy = ty;
                bodies[j].vz = tz;
            }
        }
    }
}

// function is based on previous project
uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec * (uint64_t)1000000 + tv.tv_usec;
}

void initiate_system(char* file_name) {
	FILE* fp = fopen(file_name, "r");

	fscanf(fp,"%lf%d%d", &G, &N, &T);

    bodies = (Body*)malloc(N * sizeof(Body));

	for (int i=0; i<N; i++) {
		fscanf(fp,"%lf", &bodies[i].m);
		fscanf(fp,"%lf%lf%lf", &bodies[i].px, &bodies[i].py, &bodies[i].pz);
		fscanf(fp,"%lf%lf%lf", &bodies[i].vx, &bodies[i].vy, &bodies[i].vz);
	}

	fclose(fp);
}

int main(int argc, char *argv[]) {

    if (argc != 2) {
		printf("Usage : %s <file name containing system configuration data>", argv[0]);
	} else {
		initiate_system(argv[1]);
		printf("Body   :     x              y               z           |           vx              vy              vz   ");
		for (int i=0; i<T; i++) {
			printf("\nCycle %d\n", i+1);
			compute_force();
            compute_velocity();
            compute_position();
	        // resolve_collision();
			for (int j=0; j<N; j++) {
				printf("Body %d : %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", j+1, 
				bodies[j].px, bodies[j].py, bodies[j].pz,
				bodies[j].vx, bodies[j].vy, bodies[j].vz);
			}
		}
    }
	return 0;

    // check parameter
    // if (argc != 1) {
    //     if (strcmp(argv[1], "-openmp") == 0) {
    //         printf("OpenMP approach\n");
    //     } else if (strcmp(argv[1], "-gpu") == 0) {
    //         printf("CUDA approach\n");
    //     } else {
    //         printf("parameter: -openmp | -gpu\n");
    //         exit(EXIT_FAILURE);
    //     }
    // }

    // Initialize MPI
    // MPI_Init(&argc, &argv);

    // Get MPI rank, size
    // int rank, size;
    // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // MPI_Comm_size(MPI_COMM_WORLD, &size);

    // start timing
    uint64_t start = GetTimeStamp();



    // stop timing
    uint64_t stop = GetTimeStamp();

    // master node prints result
    // if (rank == 0) {
    //     printf("Time: %llu us\n", (uint64_t) (stop - start));
    // }

    // MPI_Finalize();

    return 0;
}

/* mpicc -O3 -std=c99 -fopenmp nbody.c -o nbody */