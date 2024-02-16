/* 
 * This is just a sample code of how to use the vector instructions. 
 * This code does not get complied by the Makefile
 *
 */
void dgemm_ukernel( int m, int n, int k, double *A, int ldA,
             double *B, int ldB, double *C, int ldC )
{
  __m256d gamma_0123_0, gamma_0123_1, gamma_0123_2, gamma_0123_3;
  __m256d alpha_0123_p, beta_p_j;


  gamma_0123_0 = _mm256_loadu_pd( &gamma(0, 0) ) ;
  gamma_0123_1 = _mm256_loadu_pd( &gamma(0, 1) ) ;
  gamma_0123_2 = _mm256_loadu_pd( &gamma(0, 2) ) ;
  gamma_0123_3 = _mm256_loadu_pd( &gamma(0, 3) ) ;

  for ( int p=0; p < k; p++){
    alpha_0123_p = _mm256_loadu_pd( &alpha(0, p) ) ;

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 0) );
    gamma_0123_0 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_0 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 1) );
    gamma_0123_1 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_1 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 2) );
    gamma_0123_2 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_2 );

    beta_p_j     = _mm256_broadcast_sd( &beta( p, 3) );
    gamma_0123_3 = _mm256_fmadd_pd( alpha_0123_p, beta_p_j, gamma_0123_3 );

  }
  _mm256_storeu_pd( &gamma(0,0), gamma_0123_0 );
  _mm256_storeu_pd( &gamma(0,1), gamma_0123_1 );
  _mm256_storeu_pd( &gamma(0,2), gamma_0123_2 );
  _mm256_storeu_pd( &gamma(0,3), gamma_0123_3 );

}
