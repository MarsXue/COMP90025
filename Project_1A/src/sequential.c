#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>

#define MAX 10000
#define NOT_CONNECTED (INT_MAX)

int diameter (int distance[MAX][MAX], int nodesCount);

int distance[MAX][MAX];

/* initialize all distances to */
void Initialize(){
    for (int i=0;i<MAX;++i){
        for (int j=0;j<MAX;++j){
            distance[i][j]=NOT_CONNECTED;
        }
        distance[i][i]=0;
    }
}

uint64_t GetTimeStamp() {
    struct timeval tv;
    gettimeofday(&tv,NULL);
    return tv.tv_sec*(uint64_t)1000000+tv.tv_usec;
}

int main(){
    /* number of nodes */
    int nodeCount;

    /* Number of edges */
    int m;

    Initialize();

    /* get the node count */
    if (scanf("%d", &nodeCount) < 1) {
        fprintf (stderr, "Error reading node count\n");
        exit (1);
    }
    if (nodeCount < 1 || nodeCount > MAX) {
        fprintf (stderr, "Invalid number of nodes, %d.  Must be 1 to %d\n",
                 nodeCount, MAX);
        exit (1);
    }

    /* edge count */
    if (scanf("%d", &m) < 1) {
        fprintf (stderr, "Error reading edge count\n");
        exit (1);
    }
    if (m < nodeCount - 1 || m > nodeCount * (nodeCount - 1)) {
        fprintf (stderr, "Invalid number of edges, %d.  Must be %d to %d\n",
                 m, nodeCount-1, nodeCount * (nodeCount-1));
        exit (1);
    }

    while(m--){
        /* nodes - let the indexation begin from 0 */
        int a, b;

        /* edge weight */
        int c;

        if (scanf("%d %d %d", &a, &b, &c) < 3) {
            fprintf (stderr, "Error reading edge\n");
            exit (1);
        }
        if (a < 0 || a >= nodeCount || b < 0 || b >= nodeCount || c <= 0) {
            fprintf(stderr, "Invalid edge: from %d to %d weight %d\n", a, b, c);
            exit (1);
        }
        distance[a][b]=c;
    }

    uint64_t start = GetTimeStamp ();

    printf ("Diameter %d\n", diameter (distance, nodeCount));

    printf("Time: %ld us\n", (uint64_t) (GetTimeStamp() - start));

    return 0;
}

/******************************************************************************/
/*  Your changes here */

int diameter (int distance[MAX][MAX], int nodesCount)
{
    int i, j, k;

    for (int k=0; k<nodesCount; ++k){
        for (int i=0; i<nodesCount; ++i){
            if (distance[i][k]!=NOT_CONNECTED){
                for (int j=0; j<nodesCount; ++j){
                    if (distance[k][j]!=NOT_CONNECTED && (distance[i][j]==NOT_CONNECTED || distance[i][k]+distance[k][j]<distance[i][j])){
                        distance[i][j]=distance[i][k]+distance[k][j];
                    }
                }
            }
        }
    }

    int diameter=-1;

    /* look for the most distant pair */
    for (int i=0; i<nodesCount; ++i){
        for (int j=0; j<nodesCount; ++j){
            if (diameter<distance[i][j]){
                diameter=distance[i][j];
            }
        }
    }

    return (diameter);
}

/* The following is the exact command used to compile this code */
/* gcc -O2 graph-diameter.c -o graph-diameter */