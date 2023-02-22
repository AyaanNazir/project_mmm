#include "project2.h"


int main(int argc, char *argv[])
{
    int m, n, k;
    int m_input, n_input, k_input;
    int rsA, rsB, rsC;
    int csA, csB, csC;
    int p_begin, p_end, p_inc;
    int irep;
    int nrepeats;

  	double d_one = 1.0;
    double dtime, dtime_best;
    double diff, maxdiff = 0.0, maxreldiff=0.0, gflops;

    double *A, *B, *C, *Cold, *Cref;


    m_input = -1;
	n_input = -1;
	k_input = -1;

    nrepeats = NREPEATS;
    p_begin = P_BEGIN;
	p_end   = P_END;
	p_inc   = P_INC;

    /* Adjust first and last so that they are multiples of inc */
    printf( "%% Sweeping over matrix sizes:  %d %d %d \n", p_begin, p_end, p_inc );

    printf( "data = [\n" );
    printf( "%%   m     n     k         reference      |         current implementation \n" );
    printf( "%%                    time       GFLOPS   |    time       GFLOPS     diff \n" );
    

	for ( int p = p_begin; p <= p_end; p += p_inc ){
        maxdiff = 0.0;

        /* we will only time cases where all three matrices are square */
        if ( m_input < 0 ) m = p / ( dim_t )abs(m_input);
        else               m =     ( dim_t )    m_input;
        if ( n_input < 0 ) n = p / ( dim_t )abs(n_input);
        else               n =     ( dim_t )    n_input;
        if ( k_input < 0 ) k = p / ( dim_t )abs(k_input);
        else               k =     ( dim_t )    k_input;

        csA = m;
        csB = k;
        csC = m;

		rsA = rsB = rsC = 1;
        /* Gflops performed */
        gflops = 2.0 * m * n * k * 1e-09;

        /* Allocate space for the matrices.  We will use five arrays:
           A will be the address where A is stored.   Addressed with alpha(i,j).
           B will be the address where B is stored.   Addressed with beta(i,j).
           C will be the address where C is stored.   Addressed with gamma(i,j).

           Now, we will compute C = A B + C with via routine MyGemm
           and also with a reference implementation.  Therefore, we will
           utilize two more arrays:
 
           Cold will be the address where the original matrix C is
           stored.  

           Cref will be the address where the result of computing C = A B
           + C computed with the reference implementation will be stored.
         */

         A    = ( double * ) malloc( csA * k * sizeof( double ) );
         B    = ( double * ) malloc( csB * n * sizeof( double ) );
         C    = ( double * ) malloc( csC * n * sizeof( double ) );
         Cold = ( double * ) malloc( csC * n * sizeof( double ) );
         Cref = ( double * ) malloc( csC * n * sizeof( double ) );

         /* Generate random matrix A */
         RandomMatrix( m, k, A, rsA, csA );

         /* Generate random matrix B */
         RandomMatrix( k, n, B, rsB, csB );

         /* Generate random matrix Cold */
         RandomMatrix( m, n, Cold, rsC, csC );
    
         /* Time reference implementation provided by the BLAS library
            routine dgemm (double precision general matrix-matrix
            multiplicationn */
         for ( irep=0; irep<nrepeats; irep++ ){

             /* Copy matrix Cold to Cref */
             memcpy( Cref, Cold, csC * n * sizeof( double ) );
    
             /* start clock */
             dtime = FLA_Clock();
    
             /* Compute Cref = A B + Cref */
             bli_dgemm( BLIS_NO_TRANSPOSE, BLIS_NO_TRANSPOSE,  
					    m, n, k, &d_one, 
					    A, rsA, csA, 
					    B, rsB, csB, &d_one, 
					    Cref, rsC, csC );
             /* stop clock */
             dtime = FLA_Clock() - dtime;

             /* record the best time so far */
             if ( irep == 0 ) 
	             dtime_best = dtime;
             else
	             dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
         }
  
         printf( " %5d %5d %5d %8.4le %8.4le   ", m, n, k, dtime_best, gflops/dtime_best );
         fflush( stdout );  // We flush the output buffer because otherwise
		                    // it may throw the timings of a next
		                    // experiment.

         /* Time MyGemm */

         for ( irep=0; irep<nrepeats; irep++ ){
             /* Copy vector Cold to C */
             memcpy( C, Cold, csC * n * sizeof( double ) );
    
             /* start clock */
             dtime = FLA_Clock();
    
             /* Compute C = A B + C */
             MyGemm( m, n, k, A, rsA, csA, B, rsB, csB, C, rsC, csC );

             /* stop clock */
             dtime = FLA_Clock() - dtime;
    
             if ( irep == 0 ) 
                 dtime_best = dtime;
             else
                 dtime_best = ( dtime < dtime_best ? dtime : dtime_best );
         }

         diff    = MaxAbsDiff( m, n, C, rsC, csC, Cref, rsC, csC );
         maxdiff = ( diff > maxdiff ? diff : maxdiff );
         
         printf( " %8.4le %8.4le %8.4le \n", dtime_best, gflops/dtime_best, maxdiff );
         fflush( stdout );  // We flush the output buffer because otherwise
		                    // it may throw the timings of a next
		                    // experiment.

         /* Free the buffers */
         free( A );
         free( B );
         free( C );
         free( Cold );
         free( Cref );

    }
    printf( "];\n\n" );
    printf( "%% Maximum difference between reference and your implementation: %le.\n", maxdiff );
  
    exit( 0 );
}
