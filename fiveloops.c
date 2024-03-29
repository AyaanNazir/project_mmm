#include "project.h"

void packPortionA(int m, int k, double *A, double *packA, int csA, int rsA) {
  for (int p = 0; p < k; p++) {
    for (int i = 0; i < m; i++) {
      *packA++ = alpha(i, p);
    }
    for (int i = m; i < MR; i++) {
      *packA++ = 0;
    }
  }
}

void packageA(int m, int k, double *A, double *packA, int csA, int rsA) {
  for (int i = 0; i < m; i += MR) {
    int mb = MR > m - i ? m - i : MR;
    packPortionA(mb, k, &alpha(i, 0), packA, csA, rsA);
    packA += mb * k;
  }
}

void packPortionB(int k, int n, double *B, double *packB, int csB, int rsB) {
  for (int p = 0; p < k; p++) {
    for (int j = 0; j < n; j++) {
      *packB++ = beta(p, j);
    }
    for (int j = n; j < NR; j++) { 
      *packB++ = 0;
    }
  }
}
  
void packageB(int k, int n, double *B, double *packB, int csB, int rsB) {
  for (int j = 0; j < n; j += NR) {
    int jb = NR > n - j ? n - j : NR;
    packPortionB(k, jb, &beta(0, j), packB, csB, rsB);
    packB += k * jb;
  }
}

void gemm(int k, double *A, double *B, double *C, int rsA, int rsB, int rsC, int csA, int csB, int csC)
{
  __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3, gamma_0123_4, gamma_0123_5;
  __m256d gamma_4567_0, gamma_4567_1, gamma_4567_2, gamma_4567_3, gamma_4567_4, gamma_4567_5;
  __m256d alpha_0123_p, alpha_4567_p, beta_p_j;


  // loads for 8 x 6
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
   	
  // Unravelling loop by 4s.
  for ( int p=0; p < k; p += 4){

    // loads alpha(8)
    alpha_0123_p = _mm256_loadu_pd(&A[MR * p]);
    alpha_4567_p = _mm256_loadu_pd(&A[MR * p + 4]);

    // computes for beta(6)
    beta_p_j     = _mm256_broadcast_sd(&B[NR * p]);
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
    gamma_4567_0 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_0 );

    beta_p_j     = _mm256_broadcast_sd(&B[NR * p + 1]);
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
    gamma_4567_1 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_1 );

    beta_p_j     = _mm256_broadcast_sd(&B[NR * p + 2]);
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
    gamma_4567_2 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_2 );

    beta_p_j     = _mm256_broadcast_sd(&B[NR * p + 3]);
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
    gamma_4567_3 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_3 );

    beta_p_j     = _mm256_broadcast_sd(&B[NR * p + 4]);
    gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
    gamma_4567_4 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_4 );

    beta_p_j     = _mm256_broadcast_sd(&B[NR * p + 5]);
    gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );
    gamma_4567_5 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_5 );

    if (p + 1 < k) {
      // loads alpha(8)
      alpha_0123_p = _mm256_loadu_pd(&A[MR * (p + 1)]);
      alpha_4567_p = _mm256_loadu_pd(&A[MR * (p + 1) + 4]);

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 1)]);
      gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
      gamma_4567_0 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_0 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 1) + 1]);
      gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
      gamma_4567_1 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_1 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 1) + 2]);
      gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
      gamma_4567_2 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_2 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 1) + 3]);
      gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
      gamma_4567_3 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_3 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 1) + 4]);
      gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
      gamma_4567_4 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_4 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 1) + 5]);
      gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );
      gamma_4567_5 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_5 );
    }

    if (p + 2 < k) {
      // loads alpha(8)
      alpha_0123_p = _mm256_loadu_pd(&A[MR * (p + 2)]);
      alpha_4567_p = _mm256_loadu_pd(&A[MR * (p + 2) + 4]);

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 2)]);
      gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
      gamma_4567_0 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_0 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 2) + 1]);
      gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
      gamma_4567_1 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_1 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 2) + 2]);
      gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
      gamma_4567_2 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_2 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 2) + 3]);
      gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
      gamma_4567_3 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_3 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 2) + 4]);
      gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
      gamma_4567_4 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_4 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 2) + 5]);
      gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );
      gamma_4567_5 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_5 );
    }

    if (p + 3 < k) {
      // loads alpha(8)
      alpha_0123_p = _mm256_loadu_pd(&A[MR * (p + 3)]);
      alpha_4567_p = _mm256_loadu_pd(&A[MR * (p + 3) + 4]);

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 3)]);
      gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );
      gamma_4567_0 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_0 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 3) + 1]);
      gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );
      gamma_4567_1 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_1 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 3) + 2]);
      gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );
      gamma_4567_2 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_2 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 3) + 3]);
      gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );
      gamma_4567_3 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_3 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 3) + 4]);
      gamma_0123_4 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_4 );
      gamma_4567_4 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_4 );

      beta_p_j     = _mm256_broadcast_sd(&B[NR * (p + 3) + 5]);
      gamma_0123_5 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_5 );
      gamma_4567_5 = _mm256_fmadd_pd( alpha_4567_p, beta_p_j, gamma_4567_5 );
    }
  }

  // stores(8 x 6)
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
    int mb = MR > m - i ? m - i : MR;
     // packs C if not divisible by MR or NR
     if (mb != MR || n != NR) {

      double *packC = (double *) malloc(MR * NR * sizeof(double));
      for (int j = 0; j < n; j++) {
        for (int l = 0; l < mb; l++) {
          packC[l + (j * MR)] = gamma(i + l, j);
        }
        for (int l = mb; l < MR; l++) {
          packC[l + (j * MR)] = 0;
        }
      }
      for (int j = n; j < NR; j++) {
        for (int l = 0; l < MR; l++) {
          packC[l + (j * MR)] = 0;
        }
      }
      gemm(k, &A[i * k], B, packC, rsA, rsB, 1, csA, csB, MR);
      for (int j = 0; j < n; j++) {
        for (int l = 0; l < mb; l++) {
          gamma(i + l, j) = packC[l + (j * MR)];
        }
      }
      free(packC);
    } else {
      gemm(k, &A[i * k], B, &gamma( i,0 ), rsA, rsB, rsC, csA, csB, csC);
    }
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
  double *packA = (double *) malloc(MC * KC * sizeof(double));
  for ( int i=0; i<m; i+=MC ) {
    int ib = MC > m - i ? m - i : MC; 
    // packs A for cache
    packageA(ib, k, &alpha(i, 0), packA, csA, rsA);
    LoopTwo(ib, n, k, packA, B, &gamma( i,0 ), rsA, csA, rsB, csB, rsC, csC);
  }
  free(packA);
}

void LoopFour(int m, int n, int k, double *A, double *B, double *C,
 int rsA, int csA, int rsB, int csB, int rsC, int csC)
{
  double *packB = (double *) malloc(NC * KC * sizeof(double));
  for ( int p=0; p<k; p+=KC ) {
    int pb = KC > k - p ? k - p : KC;   
    // packs B for cache
    packageB(pb, n, &beta(p, 0), packB, csB, rsB);
    LoopThree(m, n, pb, &alpha(0, p), packB, C, rsA, csA, rsB, csB, rsC, csC);
  } 
  free(packB);
}

void fiveloops( int m, int n, int k, double *A, int rsA, int csA, 
	     double *B, int rsB, int csB,  double *C, int rsC, int csC )
{
  for ( int j=0; j<n; j+=NC ) {
    int jb = NC > n - j ? n - j : NC;

    LoopFour(m, jb, k, A, &beta( 0,j ), &gamma( 0,j ), rsA, csA, rsB, csB, rsC, csC);
  } 
}