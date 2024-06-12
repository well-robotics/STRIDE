/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:40 GMT-05:00
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
  double t2465;
  double t2457;
  double t2458;
  double t2466;
  double t2470;
  double t2459;
  double t2467;
  double t2468;
  double t2446;
  double t2471;
  double t2472;
  double t2473;
  double t2507;
  double t2509;
  double t2511;
  double t2512;
  double t2513;
  double t2528;
  double t2529;
  double t2530;
  double t2469;
  double t2490;
  double t2538;
  double t2539;
  double t2540;
  double t2515;
  double t2524;
  double t2525;
  double t2526;
  double t2527;
  double t2531;
  double t2532;
  double t2533;
  double t2534;
  double t2541;
  double t2542;
  double t2543;
  double t2544;
  double t2545;
  double t2546;
  t2465 = Cos(var1[5]);
  t2457 = Cos(var1[6]);
  t2458 = Sin(var1[5]);
  t2466 = Sin(var1[6]);
  t2470 = Cos(var1[2]);
  t2459 = -1.*t2457*t2458;
  t2467 = -1.*t2465*t2466;
  t2468 = t2459 + t2467;
  t2446 = Sin(var1[2]);
  t2471 = t2465*t2457;
  t2472 = -1.*t2458*t2466;
  t2473 = t2471 + t2472;
  t2507 = -1.*t2457;
  t2509 = 1. + t2507;
  t2511 = 0.5*t2509;
  t2512 = 0.671885*t2457;
  t2513 = t2511 + t2512;
  t2528 = t2465*t2513;
  t2529 = -0.171885*t2458*t2466;
  t2530 = t2528 + t2529;
  t2469 = -1.*t2446*t2468;
  t2490 = t2470*t2468;
  t2538 = -1.*t2465*t2457;
  t2539 = t2458*t2466;
  t2540 = t2538 + t2539;
  t2515 = -0.171885*t2465*t2466;
  t2524 = t2513*t2458;
  t2525 = 0.171885*t2465*t2466;
  t2526 = t2524 + t2525;
  t2527 = t2526*t2473;
  t2531 = t2468*t2530;
  t2532 = t2457*t2458;
  t2533 = t2465*t2466;
  t2534 = t2532 + t2533;
  t2541 = t2470*t2540;
  t2542 = t2469 + t2541;
  t2543 = 0.0732367608*var2[1]*t2542;
  t2544 = t2446*t2540;
  t2545 = t2490 + t2544;
  t2546 = 0.0732367608*var2[0]*t2545;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.0732367608*(-1.*t2446*t2473 + t2490)*var2[0] + 0.0732367608*(t2469 - 1.*t2470*t2473)*var2[1])*var2[6];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(t2543 + t2546 + 0.0732367608*(t2473*(-1.*t2458*t2513 + t2515) + t2527 + t2531 + t2530*t2534)*var2[2])*var2[6];
  p_output1[6]=(t2543 + t2546 + 0.0732367608*(t2473*(-0.171885*t2457*t2458 + t2515) + t2527 + t2531 + (0.171885*t2457*t2465 + t2529)*t2534)*var2[2] + 0.0732367608*(0.171885*t2457*t2466 - 1.*t2466*t2513)*var2[5])*var2[6];
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

#include "Ce3_vec7_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
