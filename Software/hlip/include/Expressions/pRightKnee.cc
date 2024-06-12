/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:09:56 GMT-05:00
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
  double t2557;
  double t2576;
  double t2591;
  double t2615;
  double t2530;
  double t2668;
  double t2704;
  double t2706;
  double t2710;
  double t2719;
  double t2701;
  double t2702;
  double t2703;
  double t2745;
  double t2749;
  double t2750;
  double t2612;
  double t2622;
  double t2646;
  double t2671;
  double t2674;
  double t2677;
  double t2711;
  double t2724;
  double t2734;
  double t2753;
  double t2779;
  double t2782;
  double t2790;
  double t2838;
  double t2839;
  double t2840;
  t2557 = Cos(var1[3]);
  t2576 = -1.*t2557;
  t2591 = 1. + t2576;
  t2615 = Sin(var1[3]);
  t2530 = Cos(var1[2]);
  t2668 = Sin(var1[2]);
  t2704 = Cos(var1[4]);
  t2706 = -1.*t2704;
  t2710 = 1. + t2706;
  t2719 = Sin(var1[4]);
  t2701 = t2530*t2557;
  t2702 = t2668*t2615;
  t2703 = t2701 + t2702;
  t2745 = t2557*t2668;
  t2749 = -1.*t2530*t2615;
  t2750 = t2745 + t2749;
  t2612 = -0.0265*t2591;
  t2622 = -0.0695*t2615;
  t2646 = t2612 + t2622;
  t2671 = -0.0695*t2591;
  t2674 = 0.0265*t2615;
  t2677 = t2671 + t2674;
  t2711 = -0.0265*t2710;
  t2724 = -0.2375*t2719;
  t2734 = t2711 + t2724;
  t2753 = -0.2375*t2710;
  t2779 = 0.0265*t2719;
  t2782 = t2753 + t2779;
  t2790 = t2704*t2703;
  t2838 = -1.*t2557*t2668;
  t2839 = t2530*t2615;
  t2840 = t2838 + t2839;
  p_output1[0]=t2530*t2646 + t2668*t2677 + t2703*t2734 - 0.2375*(-1.*t2703*t2719 + t2704*t2750) + t2750*t2782 - 0.0265*(t2719*t2750 + t2790) + var1[0];
  p_output1[1]=-0.0293;
  p_output1[2]=-1.*t2646*t2668 + t2530*t2677 + t2703*t2782 + t2734*t2840 - 0.0265*(t2703*t2719 + t2704*t2840) - 0.2375*(t2790 - 1.*t2719*t2840) + var1[1];
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

#include "pRightKnee.hh"

namespace SymFunction
{

void pRightKnee_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
