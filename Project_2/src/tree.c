/*
* Created for COMP90025 Parallel and Multicore Computing - Project 2, 2019
* by Hanbin Li <hanbinl1>, Wenqing Xue <wenqingx>
*
* The project is to simulate the N-Body problem in a parallel manner,
* which is a well-known topic in physics and astronomy area.
*
* This file needs to work with "NlogN_code.c" to achieve a complexity of
* O(N*log(N)).
*/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>


// represent a single particle in the system
typedef struct body {
    double m;           // mass
    double px, py;      // position x, y
    double vx, vy;      // velocity x, y
    double fx, fy;      // force x, y
} body_t;


// represent a single node of the tree
typedef struct node {
    double totalmass;
    double centerx, centery;
    double xmin, xmax;
    double ymin, ymax;
    double diag;
    struct body *body;
    struct node *NW;
    struct node *NE;
    struct node *SW;
    struct node *SE;
} node_t;


typedef enum quadrant {
    NW, NE, SW, SE
} quadrant_t;




// create a dial (xmin, xmax, ymin, ymax) and determine which sub-team the node belongs to
quadrant_t get_quadrant(double x, double y, double xmin, double xmax, double ymin, double ymax) {
    
    double midx = xmin + 0.5 * (xmax - xmin);
    double midy = ymin + 0.5 * (ymax - ymin);
    
    if (y > midy) {
        if (x > midx) {
            return NW;
        } else {
            return NE;
        }
    } else {
        if (x > midx) {
            return SE;
        } else {
            return SW;
        }
    }
}


// update the center mass after each insertion
void update_center(node_t * node, body_t * body) {
    node->centerx = (node->totalmass * node->centerx + body->m * body->px) / (node->totalmass + body->m);
    node->centery = (node->totalmass * node->centery + body->m * body->py) / (node->totalmass + body->m);
    node->totalmass += body->m;
}


// create a leaf node to insert in the quad-tree
node_t *create_node(body_t *body, double xmin, double xmax, double ymin, double ymax) {
    node_t *root = malloc(sizeof(node_t));
    
    root->totalmass = body->m;
    root->centerx = body->px;
    root->centery = body->py;
    root->xmin = xmin;
    root->xmax = xmax;
    root->ymin = ymin;
    root->ymax = ymax;
    
    root->diag = sqrt((pow(xmax - xmin, 2) + pow(ymax - ymin, 2)));
    
    root->body = body;
    
    root->NW = NULL;
    root->NE = NULL;
    root->SW = NULL;
    root->SE = NULL;

    return root;
}


// insert a paticle in the quad-tree, transforming a leaf node in a branch
void insert_body(body_t *body, node_t *node){
    
    quadrant_t existingquad, newquad;
    
    double xmid = node->xmin + 0.5 * (node->xmax - node->xmin);
    double ymid = node->ymin + 0.5 * (node->ymax - node->ymin);
    
    // first check if there is an existing quad tree
    if (node->body) {
        existingquad = get_quadrant(node->body->px, node->body->py, node->xmin, node->xmax, node->ymin, node->ymax);
        
        // insert a particle into the tree
        switch (existingquad) {
            case NW:
                node->NW = create_node(node->body, xmid, node->xmax, ymid, node->ymax);
                break;
            case NE:
                node->NE = create_node(node->body, node->xmin, xmid, ymid, node->ymax);
                break;
            case SW:
                node->SW = create_node(node->body, node->xmin, xmid, node->ymin, ymid);
                break;
            case SE:
                node->SE = create_node(node->body, xmid, node->xmax, node->ymin, ymid);
                break;
        }
        node->body = NULL;
    }

    newquad = get_quadrant(body->px, body->py, node->xmin, node->xmax, node->ymin, node->ymax);
    
    // update the centre mass of the tree
    update_center(node, body);
    
    //insert a new point in a new quadrant if it is empty, otherwise call recursively insert_body
    switch (newquad){
        case NW:
            if (!node->NW) {
                node->NW = create_node(body, xmid, node->xmax, ymid, node->ymax);
            } else {
                insert_body(body, node->NW);
            }
            break;
        case NE:
            if (!node->NE) {
                node->NE = create_node(body, node->xmin, xmid, ymid, node->ymax);
            } else {
                insert_body(body, node->NE);
            }
            break;
        case SW:
            if (!node->SW) {
                node->SW = create_node(body, node->xmin, xmid, node->ymin, ymid);
            } else {
                insert_body(body, node->SW);
            }
            break;
        case SE:
            if (!node->SE) {
                node->SE = create_node(body, xmid, node->xmax, node->ymin, ymid);
            } else {			
                insert_body(body, node->SE);
            }
            break;
    }
}


// delete the tree after use
void delete_tree(node_t *node){
    if (!node) {
        if (!node->NW) delete_tree(node->NW);
        if (!node->NE) delete_tree(node->NE);
        if (!node->SW) delete_tree(node->SW);
        if (!node->SE) delete_tree(node->SE);

        free(node);
    }
}


// sum the forces on body p
void tree_sum(node_t *node, body_t *body, double G, double treeratio){

    double factor;
    
    double dx = node->centerx - body->px;
    double dy = node->centery - body->py;
    
    double distance = sqrt(pow(dx, 2) + pow(dy, 2));
    
    if((((distance / node->diag) > treeratio) || (node->body)) && (node->body != body)) {
        factor = G * node->totalmass * body->m / pow(distance, 3);
        
        body->fx += factor * dx;
        body->fy += factor * dy;
    } else {
        if (node->NW) tree_sum(node->NW, body, G, treeratio);
        if (node->NE) tree_sum(node->NE, body, G, treeratio);
        if (node->SW) tree_sum(node->SW, body, G, treeratio);
        if (node->SE) tree_sum(node->SE, body, G, treeratio);
    }
}
