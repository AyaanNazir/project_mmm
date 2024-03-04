#ifndef MYGEMMH_
#define MYGEMMH_

/*
  Any #define or function declaration must be provided in this header.
*/

#define MR 4
#define NR 4
#define MC 128
#define NC 128
#define KC 128

void fiveloops( int, int, int, double *, int, int, double *, int, int,  double *, int, int );

#endif
