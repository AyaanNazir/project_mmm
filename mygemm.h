#ifndef MYGEMMH_
#define MYGEMMH_

/*
  Any #define or function declaration must be provided in this header.
*/

#define MR 8
#define NR 6
#define MC 96
#define NC 16384
#define KC 256

void fiveloops( int, int, int, double *, int, int, double *, int, int,  double *, int, int );

#endif
