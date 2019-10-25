/*
* Created for COMP90025 Parallel and Multicore Computing - Project 2, 2019
* by Hanbin Li <hanbinl1>, Wenqing Xue <wenqingx>
*
* The project is to simulate the N-Body problem in a parallel manner,
* which is a well-known topic in physics and astronomy area.
*
* This algorithm is implemented using Barnes Hut algorithm with OpenMP+MPI.
* This approach has a complexity of O(nlogn).
*/



#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>

#define DELTA_TIME 0.01 // delta time
#define XBOUND 1.0e6      // Width of space
#define YBOUND 1.0e6      // Height of space
#define ZBOUND 1.0e6      // Depth of space
#define RBOUND 10         // Upper bound on radius
#define THETA 1.0

//typedef struct {
//    double m;           // mass
//    double fx, fy, fz;  // force
//    double px, py, pz;  // position
//    double vx, vy, vz;  // velocity
//} Body;


/* Positional values of particle in space */
typedef struct {

   double px, py, pz;

} Position;


/* Components of velocity of particle */
typedef struct {

   double vx, vy, vz;

} Velocity;


/* Components of force experienced by particle */
typedef struct {

   double fx, fy, fz;

} Force;



/* Cubic cell representing tree node in Barnes-Hut algo. */
typedef struct Cell  {

   int index;                    // Index into arrays to identify particle's
                                 // position and mass
   int no_subcells;              // Indicate whether cell is leaf or not
   double mass;                  // Mass of particle of total mass of subtree
   double x, y, z;               // Location of cell(cube) in space
   double cx, cy, cz;            // Location of center of mass of cell
   double width, height, depth;  // Width, Height, and Depth of cell
   struct Cell* subcells[8];     // Pointers to child nodes

} Cell;


double G;
double delta_t;
int N, T;
int rank, size;
int part_size, pindex;
//Body *bodies;
Cell *root_cell;  // Root of BH octtree

Position* position;   // Current positions for all particles
Velocity* ivelocity;  // Initial velocity for all particles
Velocity* velocity;   // Velocity of particles in current processor
double* mass;         // Mass of each particle
double* radius;       // Radius of each particle
Force* force;         // Force experienced by all particles


MPI_Datatype MPI_POSITION;
MPI_Datatype MPI_VELOCITY;


// Computes the distance between two particles
double compute_distance(Position a, Position b) {
    return sqrt(pow((a.px - b.px), 2) + pow((a.py - b.py), 2) + pow((a.pz - b.pz), 2));
}

/*
 * Checks if two particles in space intersect/collide
 *
 */
int check_collision(int index1, int index2) {

   if (pow((position[index1].px - position[index2].px), 2.0) +
       pow((position[index1].py - position[index2].py), 2.0) +
       pow((position[index1].pz - position[index2].py), 2.0) <
       pow((radius[index1] + radius[index2]), 2.0)) {
       
       // Collision detected
       return 1;

   }
   return 0;
}

void reinitialize_radius() {

   int i, j;
   
   for (i = 0; i < N; i++) {

      for (j = i + 1; j < N; j++) {

         if (check_collision(i, j)) {
            double d = compute_distance(position[i], position[j]);
            radius[i] = radius[j] = d/2.0;
         }
      }
   }
}


// Computes the force for each particle
void compute_force() {
    for (int i=0; i<part_size; i++) {
        force[i].fx = 0.0;
        force[i].fy = 0.0;
        force[i].fz = 0.0;
        for (int j=0; j<N; j++) {
            if (j == (i + pindex)) continue;
            double d = compute_distance(position[i+pindex], position[j]);
            double f = (G * mass[i+pindex] * mass[j]) / (pow(d, 2.0));
            
            force[i].fx += f * ((position[j].px - position[i+pindex].px) / d);
            force[i].fy += f * ((position[j].py - position[i+pindex].py) / d);
            force[i].fz += f * ((position[j].pz - position[i+pindex].pz) / d);
        }
    }
}

// Computes the velocity for each particle
void compute_velocity() {
    for (int i=0; i<part_size; i++) {
        velocity[i].vx += (force[i].fx / mass[i+pindex]) * DELTA_TIME;
        velocity[i].vy += (force[i].fy / mass[i+pindex]) * DELTA_TIME;
        velocity[i].vz += (force[i].fz / mass[i+pindex]) * DELTA_TIME;
    }
}

// Computes the position for each particle
void compute_position() {
    for (int i=0; i<part_size; i++) {
        position[i+pindex].px += velocity[i].vx * DELTA_TIME;
        position[i+pindex].py += velocity[i].vy * DELTA_TIME;
        position[i+pindex].pz += velocity[i].vz * DELTA_TIME;
        
        // Check if particles attempt to cross boundary
           if ((position[i + pindex].px + radius[i + pindex]) >= XBOUND ||
               (position[i + pindex].px - radius[i + pindex]) <= 0)
              velocity[i].vx *= -1;
           else if ((position[i + pindex].py + radius[i + pindex] >= YBOUND) ||
                    (position[i + pindex].py - radius[i + pindex]) <= 0)
              velocity[i].vy *= -1;
           else if ((position[i + pindex].pz + radius[i + pindex]) >= ZBOUND ||
                    (position[i + pindex].pz - radius[i + pindex]) <= 0)
              velocity[i].vz *= -1;
    
    }
}


/*
 * Creates a cell (node) for use in the octtree
 *
 */
Cell* BH_create_cell(double width, double height, double depth) {

   Cell* cell = malloc(sizeof(Cell));
   cell->mass = 0;
   cell->no_subcells = 0;
   cell->index = -1;
   cell->cx = 0;
   cell->cy = 0;
   cell->cz = 0;
   cell->width = width;
   cell->height = height;
   cell->depth = depth;
   return cell;
}



/*
 * Sets the location of the subcells relative to the current cell
 *
 */
void BH_set_location_of_subcells(Cell* cell, double width, double heigth, double depth){

   // Set location of new cells
   cell->subcells[0]->x = cell->x;
   cell->subcells[0]->y = cell->y;
   cell->subcells[0]->z = cell->z;

   cell->subcells[1]->x = cell->x + width;
   cell->subcells[1]->y = cell->y;
   cell->subcells[1]->z = cell->z;

   cell->subcells[2]->x = cell->x + width;
   cell->subcells[2]->y = cell->y;
   cell->subcells[2]->z = cell->z + depth;

   cell->subcells[3]->x = cell->x;
   cell->subcells[3]->y = cell->y;
   cell->subcells[3]->z = cell->z + depth;
    
   cell->subcells[4]->x = cell->x;
   cell->subcells[4]->y = cell->y + heigth;
   cell->subcells[4]->z = cell->z;

   cell->subcells[5]->x = cell->x + width;
   cell->subcells[5]->y = cell->y + heigth;
   cell->subcells[5]->z = cell->z;

   cell->subcells[6]->x = cell->x + width;   // Coordinates of this cell marks
   cell->subcells[6]->y = cell->y + heigth;  // the mid-point of the parent cell
   cell->subcells[6]->z = cell->z + depth;   //
   
   cell->subcells[7]->x = cell->x;
   cell->subcells[7]->y = cell->y + heigth;
   cell->subcells[7]->z = cell->z + depth;

}

/*
 * Generates new subcells for the current cell, forming
 * a subtree. The current cell will no longer be a leaf
 *
 */
void BH_generate_subcells(Cell* cell) {
   
   // Calculate subcell dimensions
   double width  = cell->width / 2.0;
   double height = cell->height / 2.0;
   double depth  = cell->depth / 2.0;

   // Cell no longer a leaf
   cell->no_subcells = 8;
   
   // Create and initialize new subcells
   int i;
   for (i = 0; i < cell->no_subcells; i++) {
      cell->subcells[i] = BH_create_cell(width, height, depth);
   }
   
   BH_set_location_of_subcells(cell, width, height, depth);
}

int BH_locate_subcell(Cell* cell, int index) {

   // Determine which subcell to add the body to
   if (position[index].px > cell->subcells[6]->x){
      if (position[index].py > cell->subcells[6]->y){
         if (position[index].pz > cell->subcells[6]->z)
            return 6;
         else
            return 5;
      }
      else{
         if (position[index].pz > cell->subcells[6]->z)
            return 2;
         else
            return 1;
      }
   }
   else{
      if (position[index].py > cell->subcells[6]->y){
         if (position[index].pz > cell->subcells[6]->z)
            return 7;
         else
            return 4;
      }
      else{
         if (position[index].pz > cell->subcells[6]->z)
            return 3;
         else
            return 0;
      }
   }
}

/*
 * Added a particle to the cell. If a particle already
 * exists, the cube/cell is sub-divided adding the existing
 * and new particle to the sub cells
 *
 */
void BH_add_to_cell(Cell* cell, int index) {

   if (cell->index == -1) {
      cell->index = index;
      return;
   }
         
   BH_generate_subcells(cell);

   // The current cell's body must now be re-added to one of its subcells
   int sc1 = BH_locate_subcell(cell, cell->index);
   cell->subcells[sc1]->index = cell->index;

   // Locate subcell for new body
   int sc2 = BH_locate_subcell(cell, index);

   if (sc1 == sc2)
      BH_add_to_cell(cell->subcells[sc1], index);
   else
      cell->subcells[sc2]->index = index;
}



/*
 * Generates the octtree for the entire system of
 * particles
 *
 */
void BH_generate_octtree() {
   
   // Initialize root of octtree
   root_cell = BH_create_cell(XBOUND, YBOUND, ZBOUND);
   root_cell->index = 0;
   root_cell->x = 0;
   root_cell->y = 0;
   root_cell->z = 0;
   
   int i;
   for (i = 1; i < N; i++) {

      Cell* cell = root_cell;

      // Find which node to add the body to
      while (cell->no_subcells != 0){
         int sc = BH_locate_subcell(cell, i);
         cell = cell->subcells[sc];
      }

      BH_add_to_cell(cell, i);
   }
}


/*
 * Computes the total mass and the center of mass of
 * the current cell
 *
 */
Cell* BH_compute_cell_properties(Cell* cell){
   
   if (cell->no_subcells == 0) {
      if (cell->index != -1){
         cell->mass = mass[cell->index];
         return cell;
      }
   }
   else {
      int i;
   
      double tx = 0, ty = 0, tz = 0;
      for (i = 0; i < cell->no_subcells; i++) {
         Cell* temp = BH_compute_cell_properties(cell->subcells[i]);
         if (temp != NULL) {
            cell->mass += temp->mass;
            tx += position[temp->index].px * temp->mass;
            ty += position[temp->index].py * temp->mass;
            tz += position[temp->index].pz * temp->mass;
         }
      }
      
      // Compute center of mass
      cell->cx = tx / cell->mass;
      cell->cy = ty / cell->mass;
      cell->cz = tz / cell->mass;
   
      return cell;
   }
   return NULL;
}


/*
 * Computes the force experienced between a particle and
 * a cell
 *
 */
void BH_compute_force_from_cell(Cell* cell, int index) {
   double d = compute_distance(position[index], position[cell->index]);

   // Compute grativational force according to Newtonian's law
   double f = (G * (mass[index] * mass[cell->index]) /
                   (pow(d, 2.0)));

   // Resolve forces in each direction
   force[index - pindex].fx += f * ((position[cell->index].px - position[index].px) / d);
   force[index - pindex].fy += f * ((position[cell->index].py - position[index].py) / d);
   force[index - pindex].fz += f * ((position[cell->index].pz - position[index].pz) / d);
}




/*
 * Computes the force between the particles in the system,
 * using the clustering-approximation for long distant forces
 *
 */
void BH_compute_force_from_octtree(Cell* cell, int index) {
   
   if (cell->no_subcells == 0) {
      if (cell->index != -1 && cell->index != index) {
         BH_compute_force_from_cell(cell, index);
      }
   }
   else {
      double d = compute_distance(position[index], position[cell->index]);
      
      if (THETA > (cell->width / d)){
         // Use approximation
         BH_compute_force_from_cell(cell, index);
      }
      else {
         int i;
         for (i = 0; i < cell->no_subcells; i++) {
            BH_compute_force_from_octtree(cell->subcells[i], index);
         }
      }
   }
}




/*
 * Computes the forces experienced by the particles in space
 * using the Barnes-Hut algorithm, resulting in NlogN running
 * time
 */
void BH_compute_force(){

   int i, j;

   for (i = 0; i < part_size; i++) {

      force[i].fx = 0.0;
      force[i].fy = 0.0;
      force[i].fz = 0.0;

      BH_compute_force_from_octtree(root_cell, i + pindex);
   }
}



/*
 * Prints spaces for formatting tree output
 */
void BH_print_spaces(int number){
   int i;
   for (i = 0; i < number; i++)
      printf("  ");
}


/*
 * Prints the octtree to the console. Used for debugging.
 */
void BH_print_octtree_ex(Cell* cell, int level, int cell_no) {

   BH_print_spaces(level);
   printf("Level = %d, subcell = %d, ", level, cell_no);

   int i;

   if (cell->no_subcells == 0 && cell->index != -1) {
      BH_print_spaces(level);
      printf("position[%d] = %.2f, %.2f, %.2f; cell-location = %.2f, %.2f, %.2f, mass = %.2f;\n",
              cell->index, position[cell->index].px, position[cell->index].py,
              position[cell->index].pz, cell->x, cell->y, cell->z, mass[cell->index]);
   }
   else {
      printf("Total mass = %.2d %.2f\n", level, cell->mass);
   }

   
   if (cell->no_subcells != 0){
      level++;
      for (i = 0; i < 8; i++) {
         BH_print_octtree_ex(cell->subcells[i], level, i);
      }
   }
}


/*
 * Prints the octtree starting at the root. Used
 * for debugging
 */
void BH_print_octtree(Cell* cell){
   BH_print_octtree_ex(cell, 0, 0);
}


/*
 * Deletes the octtree
 */
void BH_delete_octtree(Cell* cell) {
   
   if (cell->no_subcells == 0) {
      free(cell);
      return;
   }

   int i;
   for (i = 0; i < cell->no_subcells; i++) {
      BH_delete_octtree(cell->subcells[i]);
   }

   free(cell);
}


/*
 * Initializes the velocity array used by each of the processes
 */
void init_velocity(){
   int i;
   for (i = 0; i < part_size; i++){
      velocity[i].vx = 0;
      velocity[i].vy = 0;
      velocity[i].vz = 0;
   }
}


void run_simulation() {
    
    if (rank == 0){
        printf("\nRunning simulation for %d bodies with %d iterations, and DELTAT = %f..\n\n", N, T, DELTA_TIME);
    }
        
    
       // Broadcast mass and position to all members in the group
    MPI_Bcast(mass, N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(position, N, MPI_POSITION, 0, MPI_COMM_WORLD);
    MPI_Scatter(ivelocity, part_size, MPI_VELOCITY, velocity, part_size, MPI_VELOCITY, 0, MPI_COMM_WORLD);


    int i;
    for (i=0; i<T; i++) {
        BH_generate_octtree();
        
        BH_compute_cell_properties(root_cell);
        BH_compute_force();
        BH_delete_octtree(root_cell);
        
        compute_velocity();
        compute_position();
        MPI_Allgather(position + (rank * part_size), part_size, MPI_POSITION,
        position, part_size, MPI_POSITION, MPI_COMM_WORLD);
    }
    
    if (rank==0) {
        //BH_print_octtree(root_cell);
        for (i = 0; i < N; i++) {
            printf("px=%f, py=%f, pz=%f\n", position[i].px, position[i].py,
                    position[i].pz);
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

   

	for (int i=0; i<N; i++) {
		fscanf(fp,"%lf%lf", &mass[i], &radius[i]);
		fscanf(fp,"%lf%lf%lf", &position[i].px, &position[i].py, &position[i].pz);
		fscanf(fp,"%lf%lf%lf", &ivelocity[i].vx, &ivelocity[i].vy, &ivelocity[i].vz);
        printf("Body %d : %lf\t%f\t%lf\t|\t%lf\t%lf\t%lf\n", i,
                    position[i].px, position[i].py, position[i].pz,
                    ivelocity[i].vx, ivelocity[i].vy, ivelocity[i].vz);
        
	}

	fclose(fp);
}




int main(int argc, char *argv[]) {
    
    MPI_Init(&argc, &argv);
    if (argc != 6) {
		printf("Usage : %s <file name containing system configuration data>", argv[0]);
	}
    
    N = strtol(argv[1], NULL, 10);
    T = strtol(argv[2], NULL, 10);
    G = strtod(argv[3], NULL);
    delta_t = strtod(argv[4], NULL);
    
    printf("N: %d\n", N);
    printf("T: %d\n", T);
    printf("G: %.2f\n", G);
    printf("delta: %.2f\n", delta_t);
    
    
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
 
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_VELOCITY);
    MPI_Type_contiguous(3, MPI_DOUBLE, &MPI_POSITION);
    MPI_Type_commit(&MPI_POSITION);
    MPI_Type_commit(&MPI_VELOCITY);
    
    
    
    part_size = N / size;

    
    pindex = rank * part_size;

    
    // Allocate memory for mass, disance, velocity and force arrays
    mass = (double *) malloc(N * sizeof(double));
    radius = (double *) malloc(N * sizeof(double));
    position = (Position *) malloc(N * sizeof(Position));
    ivelocity = (Velocity *) malloc(N * sizeof(Velocity));
    velocity = (Velocity *) malloc(part_size * sizeof(Velocity));
    force = (Force *) malloc(part_size * sizeof(Force));
    
    

    // Initialize velocity array for each process
    init_velocity();
    
    if (rank==0) {

        initiate_system(argv[5]);


    }

    
    
    
    run_simulation();
    
    
    MPI_Finalize();

    
	return 0;

   

}

/* mpicc -O3 -std=c99 -fopenmp nbody.c -o nbody */
