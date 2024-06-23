/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:55 GMT-05:00
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
  double t2444;
  double t2479;
  double t2487;
  double t2552;
  double t2412;
  double t2600;
  double t2657;
  double t2662;
  double t2667;
  double t2670;
  double t2684;
  double t2688;
  double t2696;
  double t2646;
  double t2651;
  double t2653;
  double t2530;
  double t2557;
  double t2576;
  double t2612;
  double t2615;
  double t2616;
  double t2668;
  double t2671;
  double t2674;
  double t2701;
  double t2702;
  double t2703;
  double t2705;
  double t2754;
  double t2766;
  double t2767;
  t2444 = Cos(var1[5]);
  t2479 = -1.*t2444;
  t2487 = 1. + t2479;
  t2552 = Sin(var1[5]);
  t2412 = Sin(var1[2]);
  t2600 = Cos(var1[2]);
  t2657 = Cos(var1[6]);
  t2662 = -1.*t2657;
  t2667 = 1. + t2662;
  t2670 = Sin(var1[6]);
  t2684 = t2600*t2444;
  t2688 = -1.*t2412*t2552;
  t2696 = t2684 + t2688;
  t2646 = t2444*t2412;
  t2651 = t2600*t2552;
  t2653 = t2646 + t2651;
  t2530 = -0.0695*t2487;
  t2557 = -0.0265*t2552;
  t2576 = t2530 + t2557;
  t2612 = -0.0265*t2487;
  t2615 = 0.0695*t2552;
  t2616 = t2612 + t2615;
  t2668 = -0.2375*t2667;
  t2671 = -0.0265*t2670;
  t2674 = t2668 + t2671;
  t2701 = -0.0265*t2667;
  t2702 = 0.2375*t2670;
  t2703 = t2701 + t2702;
  t2705 = t2657*t2696;
  t2754 = -1.*t2444*t2412;
  t2766 = -1.*t2600*t2552;
  t2767 = t2754 + t2766;
  p_output1[0]=t2412*t2576 + t2600*t2616 + t2653*t2674 - 0.2375*(t2653*t2657 + t2670*t2696) + t2696*t2703 - 0.0265*(-1.*t2653*t2670 + t2705) + var1[0];
  p_output1[1]=0.0293;
  p_output1[2]=t2576*t2600 - 1.*t2412*t2616 + t2674*t2696 + t2703*t2767 - 0.0265*(-1.*t2670*t2696 + t2657*t2767) - 0.2375*(t2705 + t2670*t2767) + var1[1];
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
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var1 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 1, (mwSize) 3, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "pLeftKnee.hh"

namespace SymFunction
{

void pLeftKnee_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
