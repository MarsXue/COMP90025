//
// randomgraph.c
//
// Generate a random graph on stdout.
//
// Graphs are generated using Gilbert variant of ErdÅ‘sâ€“RÃ©nyi model.
// That is, each each possible edge is present with fixed probability,
// independently.  Edges that are present are assigned an integer weight
// (distance) uniformly distributed over a fixed range from 1 up to
// some maximum weight.
//
// The resulting graph is not guaranteed to be connected, but except
// for very low edge probabilities will almost certainly be connected.
//
// Originally written by Les Kitchen <comp90025@po.ljk.id.au> for
// the subject COMP90025 Parallel and Multicore Computing
// Semester 2, 2019.
// Copyright rests with the University of Melbourne, as work-for-hire.
// Students may use this code for their study in COMP90025.
//

// Aside from using // comments, this is pretty straight C.

// The original version used drand48_r() and friends, for better
// quality of random numbers, and re-entrancy for potential
// future multi-threaded version.  But this seems not to work on
// MacOS nor Windows, nor on Spartan (which uses a much older
// gcc).  So I'm #ifdeffing code that uses older rand() interface
// as an alternative.  Random numbers won't be of as good quality,
// but should be good enough.  Hey, we're not after a high-quality
// statistical simulation; we just want to generate some random-enough
// graphs.

// Define RAND48_R to get former behavior (using rand48-style),
// probably best via -DRAND48_R on compilation command line.
// Otherwise behavior will fall back on older rand-style.
// Note that, because this will generate a different sequence of
// random numbers, the program will generate different graphs
// for the same commandline arguments, even the same seed,
// depending on compilation flags. This is somewhat unsatisfactory.
// The problem can be somewhat ameliorated by giving the executable
// program different names depending on how it was compiled, maybe
// 'randomgraph' versus 'randomgraph48_r'.

// For older versions of gcc (as on Spartan) you will need the
// -std=c99 option, to allow for-loop initial declarations.
// For very old versions of gcc, you're on your own.  This is
// the 21st century.

// Using enums instead of #define, for better syntactic behavior.

// Some weak polymorphism, to allow some flexibility in choosing
// matrix content.

#include <stdio.h>
#include <sys/time.h>
#include <stdint.h>
#include <stdlib.h>
#include <limits.h>

// Type of matrix element for this program.
typedef int sm_val_t;
#include "squarematrix.h"

enum { NOT_CONNECTED = INT_MAX };

void
usage( char *prog, size_t status )
{
#ifdef RAND48_R
  fprintf(stderr, "Usage: %s <size:size_t> <probedge:double> <maxwgt:int> <seed:long int>\n", prog);
#else
  fprintf(stderr, "Usage: %s <size:size_t> <probedge:double> <maxwgt:int> <seed:unsigned int>\n", prog);
#endif
  exit(status);
}


int main( int argc, char **argv )
{
  // Generated matrix and its size.
  sm_squarematrix_t m;
  size_t n;
  // Probability of an edge.
  double probedge;
#ifdef RAND48_R
  // Start seed for srand48_r(), and RNG state
  long int seed;
  struct drand48_data rng_state;
#else
  // Start seed for srand()
  unsigned int seed;
#endif
  // Range of weights (distances) on generated edges.
  // For this version lower is fixed and upper is argument, range is computed.
  int minwgt = 1, maxwgt, wgtrange;
#ifdef RAND48_R
  int discard;  // To swallow useless zero returned by rentrant RNG routines.
#endif
  // Count of eventually generated actual edges
  size_t edgecount = 0;

  // Four commandline args, size 'n', probability of edge 'probedge',
  // maximum weight 'maxwgt', and 'seed'.
  // Size of graph, 'n', is number of nodes.
  // For reproducibility, seed is required.  If for some
  // reason, you really want random randomness, you can
  // supply seed from some such suitable source.
  // Edges are added with probability 'probedge'.
  // This is the expected, not actual, density of edges.
  // Edges present are given an int weight uniformly distributed in range 1 to 'maxwgt'.

  // Overall argv check.
  if ( argc != 5 ) {
    usage( argv[0], 1);
  }

  // Argument: graph size 
  {
    int s = sscanf( argv[1], "%zd", &n );
    if ( s != 1 ) {
      fprintf(stderr, "%s: Bad sscanf %d on size %s giving %zd\n", argv[0], s, argv[1], n);
      usage( argv[0], 2);
    }
  }

  // Argument: probability of edge
  {
    int s = sscanf( argv[2], "%lf", &probedge );
    if ( s != 1 ) {
      fprintf(stderr, "%s: Bad sscanf %d on probedge %s giving %lf\n", argv[0], s, argv[1], probedge);
      usage( argv[0], 3);
    }
    if ( probedge < 0.0 || probedge > 1.0 ) {
      fprintf(stderr, "%s: probedge %e from %s out of range 0.0 to 1.0\n", argv[0], probedge, argv[1]);
      usage( argv[0], 4);
    }
  }

  // Argument: maxwgt, and immediately derived range
  {
    int s = sscanf( argv[3], "%d", &maxwgt );
    if ( s != 1 ) {
      fprintf(stderr, "%s: Bad sscanf %d on maxwgt %s giving %d\n", argv[0], s, argv[3], maxwgt);
      usage( argv[0], 5);
    }
    if ( maxwgt < minwgt ) {
      fprintf(stderr, "%s: maxwgt %d must be no less than minwgt %d\n", argv[0], maxwgt, minwgt);
      usage( argv[0], 5);
    }
    wgtrange = maxwgt - minwgt + 1;
  }


  // Argument: RNG seed
  {
#ifdef RAND48_R
    int s = sscanf( argv[4], "%ld", &seed );
    if ( s != 1 ) {
      fprintf(stderr, "%s: Bad sscanf %d on seed %s giving %ld\n", argv[0], s, argv[2], seed);
      usage( argv[0], 5);
    }
#else
    int s = sscanf( argv[4], "%ud", &seed );
    if ( s != 1 ) {
      fprintf(stderr, "%s: Bad sscanf %d on seed %s giving %ud\n", argv[0], s, argv[2], seed);
      usage( argv[0], 5);
    }
#endif
  }

  // Set up graph matrix, allocating storage,
  {
    sm_val_t *c = sm_alloc_content(n);
    if ( c == NULL ) {
      fprintf(stderr, "%s: Failed alloc for size %zd\n", argv[0], n);
      exit(6);
    }
    m = sm_make( n, c );
  }

#ifdef RAND48_R
  // Set up RNG.  Using re-entrant version for possible future OpenMP variant.
  discard = srand48_r( seed, &rng_state );
#else
  // Set up RNG.  Using using older version for better compatibility.
  srand( seed );
#endif

  // Randomly place weighted edges, running through all (directed) pairs of nodes.
  for ( size_t i = 0;  i < n;  i++ ) {
    for ( size_t j = 0;  j < n;  j++ ) {
      int tag;
      // With our current conventions (some implicit), a node can't be connected
      // to itself.
      if ( i == j ) {
        tag = NOT_CONNECTED;
      } else {
        double cut;
#ifdef RAND48_R
        discard = drand48_r( &rng_state, &cut );
#else
        // Note, division by floated 1+RAND_MAX to get range [0.0,1.0)
        // to be compatible with drand48
        cut = ((double) rand()) / (1.0+RAND_MAX);
#endif
        // Not sure about the < or <= here, though the difference should be negligible.
        // More a matter of aesthetics.
        // Oh, drand48 documentation says result is in [0.0,1.0).  That is, exclusive
        // at the top.  I think this is consistent with that, but I need to check.
        // I could use the excess above 'probedge' to generate the int random edge weight,
        // saving a random call.  But it could be problematic for edge probabilities insanely
        // close to 1, for which there wouldn't be enough bits left over.
        if ( cut > probedge ) {
          tag = NOT_CONNECTED;
        } else {
#ifdef RAND48_R
          discard = drand48_r( &rng_state, &cut );
#else
          // See earlier comment about 1+RAND_MAX.
          cut = ((double) rand()) / (1.0+RAND_MAX);
#endif
          tag = (int) (minwgt + cut * wgtrange);
          edgecount++;
        }
      }
      sm_set( m, i, j, tag );
    }
  }

  // Output generated graph in expected format.
  // Note, currently no checking for output errors, like
  // no-space-on-device.  We have blind faith.
  fprintf(stdout, "%zd\n%zd\n", n, edgecount);
  for ( size_t i = 0;  i < n;  i++ ) {
    for ( size_t j = 0;  j < n;  j++ ) {
      sm_val_t edge = sm_get(m,i,j);
      if ( edge != NOT_CONNECTED ) fprintf(stdout, "%zd %zd %d\n", i, j, edge);
    }
  }

  return 0;

}