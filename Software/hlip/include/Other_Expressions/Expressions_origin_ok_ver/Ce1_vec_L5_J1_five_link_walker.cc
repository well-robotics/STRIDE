/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:57 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
#include<math.h>
/**
 * Copied from Wolfram Mathematica C Definitions file mdefs.hpp
 * Changed marcos to inline functions (Eric Cousineau)
 */
inline double Power(double x, double y) { return pow(x, y); }
inline double Sqrt(double x) { return sqrt(x); }

inline double Abs(double x) { return fabs(x); }

inline double Exp(double x) { return exp(x); }
inline double Log(double x) { return log(x); }

inline double Sin(double x) { return sin(x); }
inline double Cos(double x) { return cos(x); }
inline double Tan(double x) { return tan(x); }

inline double ArcSin(double x) { return asin(x); }
inline double ArcCos(double x) { return acos(x); }
inline double ArcTan(double x) { return atan(x); }

/* update ArcTan function to use atan2 instead. */
inline double ArcTan(double x, double y) { return atan2(y,x); }

inline double Sinh(double x) { return sinh(x); }
inline double Cosh(double x) { return cosh(x); }
inline double Tanh(double x) { return tanh(x); }

const double E	= 2.71828182845904523536029;
const double Pi = 3.14159265358979323846264;
const double Degree = 0.01745329251994329576924;

inline double Sec(double x) { return 1/cos(x); }
inline double Csc(double x) { return 1/sin(x); }

#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t348;
  double t303;
  double t316;
  double t351;
  double t383;
  double t333;
  double t369;
  double t379;
  double t292;
  double t388;
  double t391;
  double t397;
  double t400;
  double t419;
  double t421;
  double t422;
  double t380;
  double t402;
  double t426;
  double t427;
  double t428;
  double t408;
  double t433;
  double t434;
  double t435;
  double t436;
  double t438;
  double t439;
  double t442;
  double t443;
  double t424;
  double t425;
  double t409;
  double t412;
  double t457;
  double t485;
  double t486;
  double t487;
  double t488;
  double t490;
  double t494;
  double t495;
  double t504;
  double t505;
  double t506;
  double t508;
  double t509;
  double t510;
  double t511;
  double t512;
  double t513;
  double t517;
  double t518;
  double t519;
  double t507;
  double t514;
  double t515;
  double t523;
  double t524;
  double t525;
  double t516;
  double t520;
  double t521;
  double t532;
  double t533;
  double t534;
  double t529;
  double t530;
  double t536;
  double t548;
  double t549;
  double t550;
  double t552;
  double t553;
  double t539;
  double t541;
  double t557;
  double t558;
  double t559;
  double t551;
  double t554;
  double t555;
  double t543;
  double t571;
  double t572;
  double t573;
  double t575;
  double t576;
  double t577;
  double t581;
  double t582;
  double t603;
  double t604;
  double t605;
  t348 = Cos(var1[5]);
  t303 = Cos(var1[6]);
  t316 = Sin(var1[5]);
  t351 = Sin(var1[6]);
  t383 = Cos(var1[2]);
  t333 = -1.*t303*t316;
  t369 = -1.*t348*t351;
  t379 = t333 + t369;
  t292 = Sin(var1[2]);
  t388 = t348*t303;
  t391 = -1.*t316*t351;
  t397 = t388 + t391;
  t400 = t383*t397;
  t419 = t303*t316;
  t421 = t348*t351;
  t422 = t419 + t421;
  t380 = t292*t379;
  t402 = t380 + t400;
  t426 = t383*t422;
  t427 = t292*t397;
  t428 = t426 + t427;
  t408 = t383*t379;
  t433 = 0.29988*t402*t428;
  t434 = -1.*t348*t303;
  t435 = t316*t351;
  t436 = t434 + t435;
  t438 = t292*t436;
  t439 = t408 + t438;
  t442 = 0.29988*t402*t439;
  t443 = t433 + t442;
  t424 = -1.*t292*t422;
  t425 = t424 + t400;
  t409 = -1.*t292*t397;
  t412 = t408 + t409;
  t457 = -1.*t292*t379;
  t485 = 0.14994*t402*t425;
  t486 = 0.14994*t412*t428;
  t487 = t383*t436;
  t488 = t457 + t487;
  t490 = 0.14994*t402*t488;
  t494 = 0.14994*t412*t439;
  t495 = t485 + t486 + t490 + t494;
  t504 = -0.022659*t303;
  t505 = -0.007367999999999986*t351;
  t506 = t504 + t505;
  t508 = -1.*t303;
  t509 = 1. + t508;
  t510 = -0.16*t509;
  t511 = -0.167368*t303;
  t512 = 0.022659*t351;
  t513 = t510 + t511 + t512;
  t517 = -1.*t316*t506;
  t518 = t348*t513;
  t519 = t517 + t518;
  t507 = t348*t506;
  t514 = t316*t513;
  t515 = t507 + t514;
  t523 = -1.*t519*t379;
  t524 = -1.*t515*t397;
  t525 = t523 + t524;
  t516 = t515*t422;
  t520 = t519*t397;
  t521 = t516 + t520;
  t532 = -1.*t348*t506;
  t533 = -1.*t316*t513;
  t534 = t532 + t533;
  t529 = 0.14994*t402*t525;
  t530 = t519*t379;
  t536 = t515*t397;
  t548 = 0.022659*t303;
  t549 = 0.007367999999999986*t351;
  t550 = t548 + t549;
  t552 = -0.007367999999999986*t303;
  t553 = t552 + t512;
  t539 = 0.14994*t521*t439;
  t541 = -1.*t515*t379;
  t557 = t348*t550;
  t558 = -1.*t316*t553;
  t559 = t557 + t558;
  t551 = t316*t550;
  t554 = t348*t553;
  t555 = t551 + t554;
  t543 = -1.*t519*t436;
  t571 = -1.*t303*t506;
  t572 = t513*t351;
  t573 = t571 + t572;
  t575 = t303*t513;
  t576 = t506*t351;
  t577 = t575 + t576;
  t581 = 0.14994*t573*t402;
  t582 = 0.14994*t577*t439;
  t603 = 0.0033974904599999994*t402;
  t604 = -0.0011047579199999977*t439;
  t605 = t603 + t604;
  p_output1[0]=var2[0]*(-0.5*(0.29988*t402*t412 + 0.29988*t425*t428)*var2[2] - 0.5*t443*var2[5] - 0.5*t443*var2[6]);
  p_output1[1]=var2[0]*(-0.5*(0.14994*Power(t412,2) + 0.14994*Power(t425,2) + 0.14994*(t409 - 1.*t383*t422)*t428 + 0.14994*t402*(-1.*t383*t397 + t457))*var2[2] - 0.5*t495*var2[5] - 0.5*t495*var2[6]);
  p_output1[2]=var2[0]*(-0.5*(0.14994*t412*t521 + 0.14994*t425*t525)*var2[2] - 0.5*(t529 + 0.14994*t402*(t422*t519 + t530 + t397*t534 + t536) + t539 + 0.14994*t428*(-1.*t397*t519 - 1.*t379*t534 + t541 + t543))*var2[5] - 0.5*(t529 + t539 + 0.14994*t428*(t541 + t543 - 1.*t397*t555 - 1.*t379*t559) + 0.14994*t402*(t530 + t536 + t422*t555 + t397*t559))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(-0.5*(0.14994*t425*t573 + 0.14994*t412*t577)*var2[2] - 0.5*(t581 + t582)*var2[5] - 0.5*(0.14994*t402*(t303*t506 - 1.*t351*t513 + t303*t550 + t351*t553) + 0.14994*t428*(t351*t550 - 1.*t303*t553 + t575 + t576) + t581 + t582)*var2[6]);
  p_output1[6]=var2[0]*(-0.5*(-0.0011047579199999977*t412 + 0.0033974904599999994*t425)*var2[2] - 0.5*t605*var2[5] - 0.5*t605*var2[6]);
}



#ifdef MATLAB_MEX_FILE

#include "mex.h"
/*
 * Main function
 */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  size_t mrows, ncols;

  double *var1,*var2;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 2)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Two input(s) required (var1,var2).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "Ce1_vec_L5_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L5_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
