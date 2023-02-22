
#include "project2.h"

void MyGemm( int m, int n, int k, double *A, int rsA, int csA,
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{

  fiveloops( m, n, k, A, rsA, csA, B, rsB, csB, C, rsC, csC);
}
  
