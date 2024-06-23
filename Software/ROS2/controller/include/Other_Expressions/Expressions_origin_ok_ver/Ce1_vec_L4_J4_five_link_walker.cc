/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:46 GMT-05:00
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
  double t269;
  double t277;
  double t313;
  double t309;
  double t285;
  double t286;
  double t287;
  double t292;
  double t294;
  double t295;
  double t272;
  double t279;
  double t280;
  double t322;
  double t307;
  double t326;
  double t327;
  double t329;
  double t282;
  double t298;
  double t303;
  double t354;
  double t356;
  double t358;
  double t333;
  double t347;
  double t348;
  double t351;
  double t359;
  double t379;
  double t380;
  double t381;
  double t412;
  double t413;
  double t402;
  double t405;
  double t408;
  double t311;
  double t316;
  double t320;
  double t382;
  double t383;
  double t386;
  double t388;
  double t391;
  double t395;
  double t360;
  double t365;
  double t431;
  double t420;
  double t421;
  double t422;
  double t321;
  double t335;
  double t441;
  double t401;
  double t409;
  double t411;
  double t415;
  double t418;
  double t442;
  double t443;
  double t444;
  double t456;
  double t457;
  double t458;
  double t461;
  double t462;
  double t463;
  double t465;
  double t466;
  double t467;
  double t459;
  double t468;
  double t472;
  double t491;
  double t492;
  double t493;
  double t487;
  double t488;
  double t489;
  double t474;
  t269 = Cos(var1[4]);
  t277 = Sin(var1[4]);
  t313 = Cos(var1[3]);
  t309 = Sin(var1[3]);
  t285 = -1.*t269;
  t286 = 1. + t285;
  t287 = -0.16*t286;
  t292 = -0.167371*t269;
  t294 = 0.022663*t277;
  t295 = t287 + t292 + t294;
  t272 = -0.022663*t269;
  t279 = -0.007370999999999989*t277;
  t280 = t272 + t279;
  t322 = Cos(var1[2]);
  t307 = Sin(var1[2]);
  t326 = t313*t269;
  t327 = -1.*t309*t277;
  t329 = t326 + t327;
  t282 = -1.*t269*t280;
  t298 = t295*t277;
  t303 = t282 + t298;
  t354 = -1.*t269*t309;
  t356 = -1.*t313*t277;
  t358 = t354 + t356;
  t333 = t322*t329;
  t347 = t269*t295;
  t348 = t280*t277;
  t351 = t347 + t348;
  t359 = t322*t358;
  t379 = t307*t358;
  t380 = t379 + t333;
  t381 = 0.14994*t303*t380;
  t412 = -0.007370999999999989*t269;
  t413 = t412 + t294;
  t402 = 0.022663*t269;
  t405 = 0.007370999999999989*t277;
  t408 = t402 + t405;
  t311 = t269*t309;
  t316 = t313*t277;
  t320 = t311 + t316;
  t382 = -1.*t313*t269;
  t383 = t309*t277;
  t386 = t382 + t383;
  t388 = t307*t386;
  t391 = t359 + t388;
  t395 = 0.14994*t351*t391;
  t360 = -1.*t307*t329;
  t365 = t359 + t360;
  t431 = -1.*t307*t358;
  t420 = -1.*t269*t413;
  t421 = t408*t277;
  t422 = t347 + t420 + t348 + t421;
  t321 = -1.*t307*t320;
  t335 = t321 + t333;
  t441 = 0.14994*t303*t365;
  t401 = t269*t280;
  t409 = t269*t408;
  t411 = -1.*t295*t277;
  t415 = t413*t277;
  t418 = t401 + t409 + t411 + t415;
  t442 = t322*t386;
  t443 = t431 + t442;
  t444 = 0.14994*t351*t443;
  t456 = -1.*t309*t280;
  t457 = t313*t295;
  t458 = t456 + t457;
  t461 = -1.*t313*t280;
  t462 = -1.*t309*t295;
  t463 = t461 + t462;
  t465 = t313*t280;
  t466 = t309*t295;
  t467 = t465 + t466;
  t459 = t458*t358;
  t468 = t467*t329;
  t472 = -1.*t467*t358;
  t491 = t313*t408;
  t492 = -1.*t309*t413;
  t493 = t491 + t492;
  t487 = t309*t408;
  t488 = t313*t413;
  t489 = t487 + t488;
  t474 = -1.*t458*t386;
  p_output1[0]=var2[3]*(-0.5*(0.14994*t303*t335 + 0.14994*t351*t365)*var2[2] - 0.5*(t381 + t395)*var2[3] - 0.5*(t381 + t395 + 0.14994*t380*t418 + 0.14994*(t320*t322 + t307*t329)*t422)*var2[4]);
  p_output1[1]=var2[3]*(-0.5*(0.14994*t303*(-1.*t320*t322 + t360) + 0.14994*t351*(-1.*t322*t329 + t431))*var2[2] - 0.5*(t441 + t444)*var2[3] - 0.5*(0.14994*t365*t418 + 0.14994*t335*t422 + t441 + t444)*var2[4]);
  p_output1[2]=var2[3]*(-0.5*(0.14994*t351*(t320*t458 + t459 + t329*t463 + t468) + 0.14994*t303*(-1.*t329*t458 - 1.*t358*t463 + t472 + t474))*var2[3] - 0.5*(0.14994*t418*(t329*t458 + t320*t467) + 0.14994*t422*(-1.*t358*t458 - 1.*t329*t467) + 0.14994*t351*(t459 + t468 + t320*t489 + t329*t493) + 0.14994*t303*(t472 + t474 - 1.*t329*t489 - 1.*t358*t493))*var2[4]);
  p_output1[3]=-0.5*(0.29988*t351*t418 + 0.29988*t303*t422)*var2[3]*var2[4];
  p_output1[4]=-0.5*(-0.0011052077399999983*t418 + 0.0033980902199999994*t422)*var2[3]*var2[4];
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

#include "Ce1_vec_L4_J4_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L4_J4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
