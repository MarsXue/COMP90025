#include <stdio.h>
#include <openssl/rand.h>

int
main (int argc, char *argv[])
{
    if (argc < 3) {
        fprintf (stderr, "usage: %s C n\n", argv[0]);
        exit (1);
    }
    
    long int C = atol (argv[1]);        /* capacity */
    int n = atoi (argv[2]);             /* number of items */
    
    int max = 4 * C / n;

    unsigned long int buff [2*n];
    RAND_seed(&max, sizeof (max));
    RAND_bytes ((void*)buff, sizeof (buff));

    printf ("%ld\n%d\n", C, n);
    int i;
    for (i = 0; i < n; i++) {
        printf ("%ld %ld\n", buff[2*i] % max, buff[2*i+1] % max);
    }
}
