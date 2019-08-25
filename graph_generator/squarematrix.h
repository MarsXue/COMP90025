//
// squarematrix.h
//
// A module to provide a thin layer of abstraction around a 2D square
// matrix, for convenience in programming things like adjacency
// matrices and edge-weight matrices for graph computations.
//
// Originally written by Les Kitchen <comp90025@po.ljk.id.au> for
// the subject COMP90025 Parallel and Multicore Computing
// Semester 2, 2019.
// Copyright rests with the University of Melbourne, as work-for-hire.
// Students may use this code for their study in COMP90025.
//

// Using enums instead of #define, for better syntactic behavior.
// Not relevant here, but left for now as message to myself.

// Some weak polymorphism, to allow some flexibility in choosing
// matrix content and indexing.

// Users of this module will need to define the type of matrix
// elements before including it, by something like:
//
//    typedef int sm_val_t;
//

#ifndef SQUAREMATRIX

#include <stdlib.h>

// Fix the index type, for which size_t will almost always serve.
// In fact, C being what it is, maybe always advertising as size_t
// is a good idea.  So, that's what I'll do.
typedef size_t sm_idx_t;

typedef struct {
  sm_idx_t sm_size;
  sm_val_t *sm_content;
} sm_squarematrix_t;

static inline sm_idx_t
sm_bytes_needed_for_content( sm_idx_t size )
{
  return size * size * sizeof (sm_val_t);

}

static inline sm_val_t
*sm_alloc_content( sm_idx_t size )
{
  return malloc( sm_bytes_needed_for_content(size) );
}

static inline void
sm_free_content( sm_squarematrix_t m)
{
  free( m.sm_content );
  return;
}

// Done this way so you have the option of having the content
// in automatic or static storage (that is, not heap allocated),
// so long as it's big enough (user's responsibility).
static inline sm_squarematrix_t
sm_make( sm_idx_t size, sm_val_t *content )
{
  sm_squarematrix_t m;
  m.sm_size = size;
  m.sm_content = content;
  return m;
}

static inline sm_idx_t
sm_size( sm_squarematrix_t m )
{
  return m.sm_size;
}

static inline sm_val_t
*sm_content( sm_squarematrix_t m)
{
  return m.sm_content;
}

static inline sm_val_t
*sm_itm_ptr( sm_squarematrix_t m, sm_idx_t row, sm_idx_t col )
{
  return &m.sm_content[ m.sm_size * row + col ];
}


static inline sm_val_t
sm_get( sm_squarematrix_t m, sm_idx_t row, sm_idx_t col )
{
  return *sm_itm_ptr(m,row,col);
 
}

static inline void
sm_set( sm_squarematrix_t m, sm_idx_t row, sm_idx_t col, sm_val_t val )
{
  *sm_itm_ptr(m,row,col) = val;
  return;
}

static inline sm_val_t
sm_setval( sm_squarematrix_t m, sm_idx_t row, sm_idx_t col, sm_val_t val )
{
  return *sm_itm_ptr(m,row,col) = val;
}

#endif


  