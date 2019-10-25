#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
    
    int N = atoi(argv[1]);

    double mass, px, py, vx, vy;

    for (int i = 0; i < N; i++) {
        mass = (double)rand() / RAND_MAX;
        px = (double)rand() / RAND_MAX * 2.0 - 1.0;
        py = (double)rand() / RAND_MAX * 2.0 - 1.0;
        vx = (double)rand() / RAND_MAX * 2.0 - 1.0;
        vy = (double)rand() / RAND_MAX * 2.0 - 1.0;
        printf("%lf %lf %lf %lf %lf\n", mass, px, py, vx, vy);
    }
    
    return 0;
}