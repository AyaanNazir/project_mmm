#include "project2.h"

void RandomMatrix( int m, int n, double *A, int rsA, int csA )
/* 
   RandomMatrix overwrite A with random values.
*/
{
  int  i, j;

  for ( i=0; i<m; i++ )
    for ( j=0; j<n; j++ )
      A[ i*rsA + j*rsA ] = drand48();
}
