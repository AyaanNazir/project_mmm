#include "project.h"

void packA(int i, int j, double *A, double *p, int csA, int rsA) {
  for (int a = 0; a < i; a += MR) {
    for (int b = 0; b < j; b++) {
      for (int c = 0; c < MR; c++) {
        *p = alpha(c, b);
        p++;
      }
    }
  }
}
  
void packB(int i, int j, double *B, double *p, int csB, int rsB) {
  for (int a = 0; a < j; a += NR) {
    for (int b = 0; b < i; b++) {
      for (int c = 0; c < NR; c++) {
        *p = beta(b, c);
        p++;
      }
    }
  }
}

void gemm(int k, double *A, double *B, double *C, int rsA, int rsB, int rsC, int csA, int csB, int csC)
{
  __m256d gamma_0123_0 = _mm256_loadu_pd( &gamma( 0,0 ) );
  __m256d gamma_0123_1 = _mm256_loadu_pd( &gamma( 0,1 ) );
  __m256d gamma_0123_2 = _mm256_loadu_pd( &gamma( 0,2 ) );
  __m256d gamma_0123_3 = _mm256_loadu_pd( &gamma( 0,3 ) );

  __m256d beta_p_j;
   	
  for ( int p=0; p<k; p++ ){
    /* load alpha( 0:3, p ) */
    __m256d alpha_0123_p = _mm256_loadu_pd( &alpha(0, p) );

    /* load beta( p, 0 ); update gamma( 0:3, 0 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );

    /* load beta( p, 1 ); update gamma( 0:3, 1 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );

    /* load beta( p, 2 ); update gamma( 0:3, 2 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );

    /* load beta( p, 3 ); update gamma( 0:3, 3 ) */
    beta_p_j = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
  }
}

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  for (int a = 0; a < n; a += NC) {
    // double *Bpack = calloc(sizeof(double) * NC * KC);
    for (int b = 0; b < k; b += KC) {
      // packB(KC, n, &beta(b, 0), Bpack, csB, rsB);
      // double *Apack = calloc(sizeof(double) * MC * KC);
      for (int c = 0; c < m; c += MC) {
        // packA(MC, k, &alpha(c, 0), csA, Apack, csA, csB);
        for (int d = 0; d < NC; d += NR) {
          for (int e = 0; e < MC; e += MR) {
            gemm(k, &alpha(d, a), &beta(a, e), &gamma(d, e), rsA, rsB, rsC, csA, csB, csC);
          }
        }
      }
      // free(Apack);
    }
    // free(Bpack);
  }
}
