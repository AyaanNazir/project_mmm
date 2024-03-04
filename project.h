#ifndef PROJECTH_
#define PROJECTH_
#define _XOPEN_SOURCE 


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include<immintrin.h>

#include "blis.h"
#include "mygemm.h"



#define dabs( x ) ( (x) < 0 ? -(x) : x )

#define alpha( i,j ) A[ (i)*rsA + (j)*csA ]   // map alpha( i,j ) to array A 
#define beta( i,j )  B[ (i)*rsB + (j)*csB ]   // map beta( i,j )  to array B
#define gamma( i,j ) C[ (i)*rsC + (j)*csC ]   // map gamma( i,j ) to array C

double FLA_Clock();      // This is a routine for extracting elapsed
			 // time borrowed from the libflame library

/* MaxAbsDiff computes the maximum absolute difference over all
   corresponding elements of two matrices */
double MaxAbsDiff( int, int, double *, int, int, double *, int, int );

/* RandomMatrix overwrites a matrix with random values */
void RandomMatrix( int, int, double *, int, int );

/* My_Gemm is a common interface to all the implementations we will 
   develop so we don't have to keep rewriting this driver routine. */
void MyGemm( int, int, int, double *, int, int, double *, int, int,  double *, int, int );

#endif
