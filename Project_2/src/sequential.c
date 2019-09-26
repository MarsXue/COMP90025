/*
 * Sequential source code for N-Body problem
 * 
 * https://rosettacode.org/wiki/N-body_problem#C
 * 
 * gcc -O3 sequential.c -o sequential
 * ./sequential input.txt
 */

#include<stdlib.h>
#include<stdio.h>
#include<math.h>

typedef struct{
	double x, y, z;
} vector;

int bodies, timeSteps;
double *masses, GravConstant;
vector *positions, *velocities, *accelerations;

vector add_vectors(vector a, vector b) {
	vector c = {a.x+b.x, a.y+b.y, a.z+b.z};
	return c;
}

vector scale_vector(double b, vector a) {
	vector c = {b*a.x, b*a.y, b*a.z};
	return c;
}

vector subtract_vectors(vector a, vector b) {
	vector c = {a.x-b.x, a.y-b.y, a.z-b.z};
	return c;
}

double mod(vector a){
	return sqrt(a.x*a.x + a.y*a.y + a.z*a.z);
}

void initiate_system(char* file_name) {
	FILE* fp = fopen(file_name, "r");

	fscanf(fp,"%lf%d%d", &GravConstant, &bodies, &timeSteps);

	masses = (double*)malloc(bodies * sizeof(double));
	positions = (vector*)malloc(bodies * sizeof(vector));
	velocities = (vector*)malloc(bodies * sizeof(vector));
	accelerations = (vector*)malloc(bodies * sizeof(vector));

	for (int i=0; i<bodies; i++) {
		fscanf(fp,"%lf", &masses[i]);
		fscanf(fp,"%lf%lf%lf", &positions[i].x, &positions[i].y, &positions[i].z);
		fscanf(fp,"%lf%lf%lf", &velocities[i].x, &velocities[i].y, &velocities[i].z);
	}

	fclose(fp);
}

void resolve_collisions() {
	for (int i=0; i<bodies-1; i++) {
		for (int j=i+1; j<bodies; j++) {
			if (positions[i].x == positions[j].x && 
				positions[i].y == positions[j].y &&
				positions[i].z == positions[j].z) {
				vector temp = velocities[i];
				velocities[i] = velocities[j];
				velocities[j] = temp;
			}
		}
	}
}

void compute_accelerations() {
	for (int i=0; i<bodies; i++) {
		accelerations[i].x = 0;
		accelerations[i].y = 0;
		accelerations[i].z = 0;
		for (int j=0; j<bodies; j++) {
			if (i != j) {
				accelerations[i] = add_vectors(accelerations[i], scale_vector(GravConstant*masses[j]/pow(mod(subtract_vectors(positions[i],positions[j])),3), subtract_vectors(positions[j],positions[i])));
			}
		}
	}
}

void compute_velocities(){
	for (int i=0; i<bodies; i++) {
		velocities[i] = add_vectors(velocities[i],accelerations[i]);
	}
}

void compute_positions(){
	for (int i=0;i<bodies; i++) {
		positions[i] = add_vectors(positions[i], add_vectors(velocities[i], scale_vector(0.5,accelerations[i])));
	}
}

void simulate(){
	compute_accelerations();
	compute_positions();
	compute_velocities();
	resolve_collisions();
}

int main(int argc,char* argv[]) {
	if (argc!=2) {
		printf("Usage : %s <file name containing system configuration data>",argv[0]);
	} else {
		initiate_system(argv[1]);
		printf("Body   :     x              y               z           |           vx              vy              vz   ");
		for (int i=0; i<timeSteps; i++) {
			printf("\nCycle %d\n", i+1);
			simulate();
			for (int j=0; j<bodies; j++) {
				printf("Body %d : %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n",j+1,positions[j].x,positions[j].y,positions[j].z,velocities[j].x,velocities[j].y,velocities[j].z);
			}
		}
	}
	return 0;
}