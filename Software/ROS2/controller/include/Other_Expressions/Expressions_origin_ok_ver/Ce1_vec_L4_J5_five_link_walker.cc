/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:49 GMT-05:00
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
  double t294;
  double t279;
  double t287;
  double t298;
  double t326;
  double t272;
  double t333;
  double t335;
  double t336;
  double t359;
  double t365;
  double t368;
  double t347;
  double t369;
  double t391;
  double t397;
  double t400;
  double t401;
  double t402;
  double t405;
  double t408;
  double t409;
  double t412;
  double t415;
  double t292;
  double t303;
  double t311;
  double t378;
  double t379;
  double t423;
  double t433;
  double t434;
  double t435;
  double t436;
  double t437;
  double t445;
  double t446;
  double t447;
  double t448;
  double t449;
  double t452;
  double t453;
  double t455;
  double t457;
  double t459;
  double t460;
  double t464;
  double t468;
  double t469;
  double t470;
  double t473;
  double t475;
  double t476;
  double t465;
  double t477;
  double t488;
  double t490;
  double t491;
  double t495;
  double t496;
  double t481;
  double t500;
  double t501;
  double t502;
  double t494;
  double t497;
  double t498;
  double t483;
  t294 = Cos(var1[3]);
  t279 = Cos(var1[4]);
  t287 = Sin(var1[3]);
  t298 = Sin(var1[4]);
  t326 = Cos(var1[2]);
  t272 = Sin(var1[2]);
  t333 = t294*t279;
  t335 = -1.*t287*t298;
  t336 = t333 + t335;
  t359 = -1.*t279*t287;
  t365 = -1.*t294*t298;
  t368 = t359 + t365;
  t347 = t326*t336;
  t369 = t326*t368;
  t391 = t272*t368;
  t397 = t391 + t347;
  t400 = 0.0033980902199999994*t397;
  t401 = -1.*t294*t279;
  t402 = t287*t298;
  t405 = t401 + t402;
  t408 = t272*t405;
  t409 = t369 + t408;
  t412 = -0.0011052077399999983*t409;
  t415 = t400 + t412;
  t292 = t279*t287;
  t303 = t294*t298;
  t311 = t292 + t303;
  t378 = -1.*t272*t336;
  t379 = t369 + t378;
  t423 = -1.*t272*t368;
  t433 = 0.0033980902199999994*t379;
  t434 = t326*t405;
  t435 = t423 + t434;
  t436 = -0.0011052077399999983*t435;
  t437 = t433 + t436;
  t445 = -0.022663*t279;
  t446 = -0.007370999999999989*t298;
  t447 = t445 + t446;
  t448 = -1.*t287*t447;
  t449 = -1.*t279;
  t452 = 1. + t449;
  t453 = -0.16*t452;
  t455 = -0.167371*t279;
  t457 = 0.022663*t298;
  t459 = t453 + t455 + t457;
  t460 = t294*t459;
  t464 = t448 + t460;
  t468 = -1.*t294*t447;
  t469 = -1.*t287*t459;
  t470 = t468 + t469;
  t473 = t294*t447;
  t475 = t287*t459;
  t476 = t473 + t475;
  t465 = t464*t368;
  t477 = t476*t336;
  t488 = 0.022663*t279;
  t490 = 0.007370999999999989*t298;
  t491 = t488 + t490;
  t495 = -0.007370999999999989*t279;
  t496 = t495 + t457;
  t481 = -1.*t476*t368;
  t500 = t294*t491;
  t501 = -1.*t287*t496;
  t502 = t500 + t501;
  t494 = t287*t491;
  t497 = t294*t496;
  t498 = t494 + t497;
  t483 = -1.*t464*t405;
  p_output1[0]=var2[4]*(-0.5*(0.0033980902199999994*(-1.*t272*t311 + t347) - 0.0011052077399999983*t379)*var2[2] - 0.5*t415*var2[3] - 0.5*t415*var2[4]);
  p_output1[1]=var2[4]*(-0.5*(0.0033980902199999994*(-1.*t311*t326 + t378) - 0.0011052077399999983*(-1.*t326*t336 + t423))*var2[2] - 0.5*t437*var2[3] - 0.5*t437*var2[4]);
  p_output1[2]=var2[4]*(-0.5*(-0.0011052077399999983*(t311*t464 + t465 + t336*t470 + t477) + 0.0033980902199999994*(-1.*t336*t464 - 1.*t368*t470 + t481 + t483))*var2[3] - 0.5*(-0.0011052077399999983*(t465 + t477 + t311*t498 + t336*t502) + 0.0033980902199999994*(t481 + t483 - 1.*t336*t498 - 1.*t368*t502))*var2[4]);
  p_output1[3]=-0.5*(0.0033980902199999994*(t298*t447 + t279*t459 + t298*t491 - 1.*t279*t496) - 0.0011052077399999983*(t279*t447 - 1.*t298*t459 + t279*t491 + t298*t496))*Power(var2[4],2);
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=0;
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

#include "Ce1_vec_L4_J5_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L4_J5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
