#pragma once
////////////////////////////////////////////////////////////////////////////////
#define _USE_MATH_DEFINES
#include <iostream>
#include <vector>
#include <string>
#include <stdio.h>
#include <time.h>
using namespace std;
#include "vecmath.h"

////////////////////////////////////////////////////////////////////////////////
#define _USE_MATH_DEFINES
#define loopi(start_l,end_l,step_l) for ( unsigned i=start_l;i<end_l;i+=step_l )
#define loopi(start_l,end_l) for ( unsigned i=start_l;i<end_l;++i )
#define loopj(start_l,end_l,step_l) for ( unsigned j=start_l;j<end_l;j+=step_l )
#define loopj(start_l,end_l) for ( unsigned j=start_l;j<end_l;++j )
#define loopk(start_l,end_l,step_l) for ( unsigned k=start_l;k<end_l;k+=step_l )
#define loopk(start_l,end_l) for ( unsigned k=start_l;k<end_l;++k )
#define loopm(start_l,end_l) for ( unsigned m=start_l;m<end_l;++m )
#define loopl(start_l,end_l) for ( unsigned l=start_l;l<end_l;++l )
#define loopn(start_l,end_l) for ( unsigned n=start_l;n<end_l;++n )
//#define loop(var_l,start_l,end_l) for ( unsigned var_l=start_l;var_l<end_l;++var_l )
#define loops(a_l,start_l,end_l,step_l) for ( a_l = start_l;a_l<end_l;a_l+=step_l )

/*
#ifndef byte
#define byte unsigned char
#endif
*/
#ifndef ushort
#define ushort unsigned short
#endif

#ifndef uint
#define uint unsigned int
#endif
/*
#ifndef uchar
#define uchar char
#endif
*/
////////////////////////////////////////////////////////////////////////////////

#ifdef _cplusplus
extern "C"{
#endif // _cplusplus
double CoreCubicInterpolate(
   double y0,double y1,
   double y2,double y3,
   double mu);

#ifdef _cplusplus
}
#endif // _cplusplus
