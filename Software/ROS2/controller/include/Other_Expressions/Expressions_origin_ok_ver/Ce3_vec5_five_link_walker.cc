/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:39 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
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


#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t2298;
  double t2283;
  double t2284;
  double t2299;
  double t2303;
  double t2285;
  double t2300;
  double t2301;
  double t2259;
  double t2304;
  double t2308;
  double t2309;
  double t2379;
  double t2388;
  double t2390;
  double t2391;
  double t2392;
  double t2410;
  double t2411;
  double t2412;
  double t2302;
  double t2338;
  double t2420;
  double t2421;
  double t2422;
  double t2394;
  double t2406;
  double t2407;
  double t2408;
  double t2409;
  double t2413;
  double t2414;
  double t2415;
  double t2416;
  double t2423;
  double t2424;
  double t2425;
  double t2426;
  double t2427;
  double t2428;
  t2298 = Cos(var1[3]);
  t2283 = Cos(var1[4]);
  t2284 = Sin(var1[3]);
  t2299 = Sin(var1[4]);
  t2303 = Cos(var1[2]);
  t2285 = -1.*t2283*t2284;
  t2300 = -1.*t2298*t2299;
  t2301 = t2285 + t2300;
  t2259 = Sin(var1[2]);
  t2304 = t2298*t2283;
  t2308 = -1.*t2284*t2299;
  t2309 = t2304 + t2308;
  t2379 = -1.*t2283;
  t2388 = 1. + t2379;
  t2390 = 0.5*t2388;
  t2391 = 0.671885*t2283;
  t2392 = t2390 + t2391;
  t2410 = t2298*t2392;
  t2411 = -0.171885*t2284*t2299;
  t2412 = t2410 + t2411;
  t2302 = -1.*t2259*t2301;
  t2338 = t2303*t2301;
  t2420 = -1.*t2298*t2283;
  t2421 = t2284*t2299;
  t2422 = t2420 + t2421;
  t2394 = -0.171885*t2298*t2299;
  t2406 = t2392*t2284;
  t2407 = 0.171885*t2298*t2299;
  t2408 = t2406 + t2407;
  t2409 = t2408*t2309;
  t2413 = t2301*t2412;
  t2414 = t2283*t2284;
  t2415 = t2298*t2299;
  t2416 = t2414 + t2415;
  t2423 = t2303*t2422;
  t2424 = t2302 + t2423;
  t2425 = 0.0732367608*var2[1]*t2424;
  t2426 = t2259*t2422;
  t2427 = t2338 + t2426;
  t2428 = 0.0732367608*var2[0]*t2427;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.0732367608*(-1.*t2259*t2309 + t2338)*var2[0] + 0.0732367608*(t2302 - 1.*t2303*t2309)*var2[1])*var2[4];
  p_output1[3]=(t2425 + t2428 + 0.0732367608*(t2309*(-1.*t2284*t2392 + t2394) + t2409 + t2413 + t2412*t2416)*var2[2])*var2[4];
  p_output1[4]=(t2425 + t2428 + 0.0732367608*(t2309*(-0.171885*t2283*t2284 + t2394) + t2409 + t2413 + (0.171885*t2283*t2298 + t2411)*t2416)*var2[2] + 0.0732367608*(0.171885*t2283*t2299 - 1.*t2299*t2392)*var2[3])*var2[4];
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

#include "Ce3_vec5_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
