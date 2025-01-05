/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:14 GMT-05:00
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

/*
 * Sub functions
 */
static void output1(double *p_output1,const double *var1)
{
  double t2610;
  double t2614;
  double t2638;
  double t2669;
  double t2551;
  double t2724;
  double t2760;
  double t2762;
  double t2763;
  double t2768;
  double t2753;
  double t2758;
  double t2759;
  double t2799;
  double t2802;
  double t2806;
  double t2659;
  double t2673;
  double t2679;
  double t2727;
  double t2728;
  double t2731;
  double t2767;
  double t2776;
  double t2781;
  double t2807;
  double t2824;
  double t2836;
  double t2840;
  double t2894;
  double t2895;
  double t2896;
  t2610 = Cos(var1[3]);
  t2614 = -1.*t2610;
  t2638 = 1. + t2614;
  t2669 = Sin(var1[3]);
  t2551 = Cos(var1[2]);
  t2724 = Sin(var1[2]);
  t2760 = Cos(var1[4]);
  t2762 = -1.*t2760;
  t2763 = 1. + t2762;
  t2768 = Sin(var1[4]);
  t2753 = t2551*t2610;
  t2758 = t2724*t2669;
  t2759 = t2753 + t2758;
  t2799 = t2610*t2724;
  t2802 = -1.*t2551*t2669;
  t2806 = t2799 + t2802;
  t2659 = -0.0265*t2638;
  t2673 = -0.0695*t2669;
  t2679 = t2659 + t2673;
  t2727 = -0.0695*t2638;
  t2728 = 0.0265*t2669;
  t2731 = t2727 + t2728;
  t2767 = -0.0265*t2763;
  t2776 = -0.2375*t2768;
  t2781 = t2767 + t2776;
  t2807 = -0.2375*t2763;
  t2824 = 0.0265*t2768;
  t2836 = t2807 + t2824;
  t2840 = t2760*t2759;
  t2894 = -1.*t2610*t2724;
  t2895 = t2551*t2669;
  t2896 = t2894 + t2895;
  p_output1[0]=t2551*t2679 + t2724*t2731 + t2759*t2781 - 0.2375*(-1.*t2759*t2768 + t2760*t2806) + t2806*t2836 - 0.0265*(t2768*t2806 + t2840) + var1[0];
  p_output1[1]=-0.0293;
  p_output1[2]=-1.*t2679*t2724 + t2551*t2731 + t2759*t2836 + t2781*t2896 - 0.0265*(t2759*t2768 + t2760*t2896) - 0.2375*(t2840 - 1.*t2768*t2896) + var1[1];
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
