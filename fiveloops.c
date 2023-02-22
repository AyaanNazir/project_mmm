#define alpha( i,j ) A[ i*rsA + j*csA ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ i*rsB + j*csB ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ i*rsC + j*csC ]   // map gamma( i,j ) to array C

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  for ( int i=0; i<m; i++ )
    for ( int j=0; j<n; j++ )
      for ( int p=0; p<k; p++ )
        gamma( i,j ) += alpha( i,p ) * beta( p,j );
}
  
