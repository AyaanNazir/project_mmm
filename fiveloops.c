#include "project.h"

void packPortionA(int k, double *A, double *packA, int csA, int rsA) {
  for (int p = 0; p < k; p++) {
    for (int i = 0; i < MR; i++) {
      *packA++ = alpha(i, p);
    }
  }
}

void packageA(int m, int k, double *A, double *packA, int csA, int rsA) {
  for (int i = 0; i < m; i += MR) {
    packPortionA(k, &alpha(i, 0), packA, csA, rsA);
    packA += MR * k;
  }
}

void packPortionB(int k, double *B, double *packB, int csB, int rsB) {
  for (int p = 0; p < k; p++) {
    for (int j = 0; j < NR; j++) {
      *packB++ = beta(p, j);
    }
  }
}
  
void packageB(int k, int n, double *B, double *packB, int csB, int rsB) {
  for (int j = 0; j < n; j += NR) {
    packPortionB(k, &beta(0, j), packB, csB, rsB);
    packB += k * NR;
  }
}

void gemm(int k, double *A, double *B, double *C, int rsA, int rsB, int rsC, int csA, int csB, int csC)
{
  __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3, gamma_0123_4, gamma_0123_5;
  __m256d gamma_4567_0, gamma_4567_1, gamma_4567_2, gamma_4567_3, gamma_4567_4, gamma_4567_5;
  __m256d alpha_0123_p, alpha_4567_p, beta_p_j;


  gamma_0123_0 = _mm256_loadu_pd( &gamma(0, 0) ) ;
  gamma_0123_1 = _mm256_loadu_pd( &gamma(0, 1) ) ;
  gamma_0123_2 = _mm256_loadu_pd( &gamma(0, 2) ) ;
  gamma_0123_3 = _mm256_loadu_pd( &gamma(0, 3) ) ;
  gamma_0123_4 = _mm256_loadu_pd( &gamma(0, 4) ) ;
  gamma_0123_5 = _mm256_loadu_pd( &gamma(0, 5) ) ;
  gamma_4567_0 = _mm256_loadu_pd( &gamma(4, 0) ) ;
  gamma_4567_1 = _mm256_loadu_pd( &gamma(4, 1) ) ;
  gamma_4567_2 = _mm256_loadu_pd( &gamma(4, 2) ) ;
  gamma_4567_3 = _mm256_loadu_pd( &gamma(4, 3) ) ;
  gamma_4567_4 = _mm256_loadu_pd( &gamma(4, 4) ) ;
  gamma_4567_5 = _mm256_loadu_pd( &gamma(4, 5) ) ;
   	
  for ( int p=0; p < k; p++){
    alpha_0123_p = _mm256_loadu_pd(A);
    alpha_4567_p = _mm256_loadu_pd(A + 4);

    beta_p_j     = _mm256_broadcast_sd(B);
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    gamma_4567_0 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_0 );

    beta_p_j     = _mm256_broadcast_sd(B + 1);
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    gamma_4567_1 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_1 );

    beta_p_j     = _mm256_broadcast_sd(B + 2);
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    gamma_4567_2 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_2 );

    beta_p_j     = _mm256_broadcast_sd(B + 3);
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
    gamma_4567_3 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_3 );

    beta_p_j     = _mm256_broadcast_sd(B + 4);
    gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
    gamma_4567_4 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_4 );

    beta_p_j     = _mm256_broadcast_sd(B + 5);
    gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );
    gamma_4567_5 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_5 );

    A += MR;
    B += NR;
  }
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );
  _mm256_storeu_pd( &gamma(0,4), gamma_0123_4 );
  _mm256_storeu_pd( &gamma(0,5), gamma_0123_5 );
  _mm256_storeu_pd( &gamma(4,0), gamma_4567_0 );
  _mm256_storeu_pd( &gamma(4,1), gamma_4567_1 );
  _mm256_storeu_pd( &gamma(4,2), gamma_4567_2 );
  _mm256_storeu_pd( &gamma(4,3), gamma_4567_3 );
  _mm256_storeu_pd( &gamma(4,4), gamma_4567_4 );
  _mm256_storeu_pd( &gamma(4,5), gamma_4567_5 );
}

void LoopOne(int m, int n, int k, double *A, double *B, double *C,
 int rsA, int csA, int rsB, int csB, int rsC, int csC)
{
  for ( int i=0; i<m; i+=MR ) {
    int ib = MR > m - i ? m - i : MR;

    gemm(k, &A[i * k], B, &gamma( i,0 ), rsA, rsB, rsC, csA, csB, csC);
  }
}

void LoopTwo(int m, int n, int k, double *A, double *B, double *C,
 int rsA, int csA, int rsB, int csB, int rsC, int csC)
{
  for ( int j=0; j<n; j+=NR ) {
    int jb = NR > n - j ? n - j : NR;

    LoopOne(m, jb, k, A, &B[j * k], &gamma( 0,j ), rsA, csA, rsB, csB, rsC, csC);
  }
}

void LoopThree(int m, int n, int k, double *A, double *B, double *C,
 int rsA, int csA, int rsB, int csB, int rsC, int csC)
{
  double *packA = (double *) calloc(MC * KC, sizeof(double));
  for ( int i=0; i<m; i+=MC ) {
    int ib = MC > m - i ? m - i : MC; 
    packageA(ib, k, &alpha(i, 0), packA, csA, rsA);
    LoopTwo(ib, n, k, packA, B, &gamma( i,0 ), rsA, csA, rsB, csB, rsC, csC);
  }
  free(packA);
}

void LoopFour(int m, int n, int k, double *A, double *B, double *C,
 int rsA, int csA, int rsB, int csB, int rsC, int csC)
{
  double *packB = (double *) calloc(NC * KC, sizeof(double));
  for ( int p=0; p<k; p+=KC ) {
    int pb = KC > k - p ? k - p : KC;   
    packageB(pb, n, &beta(p, 0), packB, csB, rsB);
    LoopThree(m, n, pb, &alpha(0, p), packB, C, rsA, csA, rsB, csB, rsC, csC);
  } 
  free(packB);
}

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  for ( int j=0; j<n; j+=NC ) {
    int jb = NC > n - j ? n - j : NC;    /* Last loop may not involve a full block */

    LoopFour(m, jb, k, A, &beta( 0,j ), &gamma( 0,j ), rsA, csA, rsB, csB, rsC, csC);
  } 
}