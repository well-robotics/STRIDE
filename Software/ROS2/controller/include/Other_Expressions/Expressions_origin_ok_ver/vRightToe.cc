/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:32 GMT-05:00
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
static void output1(double *p_output1,const double *var1)
{
  double t220;
  double t261;
  double t208;
  double t290;
  double t321;
  double t324;
  double t326;
  double t335;
  double t346;
  double t347;
  double t348;
  double t315;
  double t318;
  double t319;
  double t339;
  double t373;
  double t374;
  double t375;
  double t359;
  double t364;
  double t263;
  double t395;
  double t396;
  double t327;
  double t341;
  double t349;
  double t350;
  double t354;
  double t406;
  double t407;
  double t408;
  double t399;
  double t400;
  double t402;
  double t403;
  double t404;
  double t356;
  double t416;
  double t241;
  double t267;
  double t293;
  double t295;
  double t313;
  double t412;
  double t385;
  double t388;
  double t370;
  double t371;
  double t376;
  double t378;
  double t379;
  double t413;
  double t414;
  t220 = Cos(var1[3]);
  t261 = Sin(var1[3]);
  t208 = Sin(var1[2]);
  t290 = Cos(var1[2]);
  t321 = Cos(var1[4]);
  t324 = -1.*t321;
  t326 = 1. + t324;
  t335 = Sin(var1[4]);
  t346 = -1.*t290*t220;
  t347 = -1.*t208*t261;
  t348 = t346 + t347;
  t315 = t220*t208;
  t318 = -1.*t290*t261;
  t319 = t315 + t318;
  t339 = -0.2375*t335;
  t373 = t290*t220;
  t374 = t208*t261;
  t375 = t373 + t374;
  t359 = -1.*t319*t335;
  t364 = t321*t319;
  t263 = -0.0695*t261;
  t395 = -1.*t220;
  t396 = 1. + t395;
  t327 = -0.0265*t326;
  t341 = t327 + t339;
  t349 = -0.2375*t326;
  t350 = 0.0265*t335;
  t354 = t349 + t350;
  t406 = -1.*t220*t208;
  t407 = t290*t261;
  t408 = t406 + t407;
  t399 = -0.0265*t396;
  t400 = t399 + t263;
  t402 = -0.0695*t396;
  t403 = 0.0265*t261;
  t404 = t402 + t403;
  t356 = t321*t348;
  t416 = t321*t408;
  t241 = 0.0265*t220;
  t267 = t241 + t263;
  t293 = -0.0695*t220;
  t295 = -0.0265*t261;
  t313 = t293 + t295;
  t412 = t321*t375;
  t385 = -1.*t375*t335;
  t388 = t364 + t385;
  t370 = 0.0265*t321;
  t371 = t370 + t339;
  t376 = -0.2375*t321;
  t378 = -0.0265*t335;
  t379 = t376 + t378;
  t413 = -1.*t408*t335;
  t414 = t412 + t413;
  p_output1[0]=var1[7] + (t354*t375 - 1.*t208*t400 + t290*t404 + t341*t408 + 0.0225*t414 - 0.0265*(t335*t375 + t416))*var1[9] + (t208*t267 + t290*t313 + t319*t341 + t348*t354 + 0.0225*(t356 + t359) - 0.0265*(t335*t348 + t364))*var1[10] + (t319*t371 + 0.0225*(t359 - 1.*t321*t375) + t375*t379 - 0.0265*t388)*var1[11];
  p_output1[1]=var1[8] + (t341*t348 - 1.*t290*t400 - 1.*t208*t404 + t354*t408 - 0.0265*(t356 + t335*t408) + 0.0225*(-1.*t335*t348 + t416))*var1[9] + (t267*t290 - 1.*t208*t313 + t319*t354 + t341*t375 + 0.0225*t388 - 0.0265*(t319*t335 + t412))*var1[10] + (t371*t375 + t379*t408 + 0.0225*(t385 - 1.*t321*t408) - 0.0265*t414)*var1[11];
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

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
    }
  else if( nlhs > 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:maxlhs", "Too many output arguments.");
    }

  /*  The input must be a noncomplex double vector or scaler.  */
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  if( !mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) ||
    ( !(mrows == 14 && ncols == 1) && 
      !(mrows == 1 && ncols == 14))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "vRightToe.hh"

namespace SymFunction
{

void vRightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
