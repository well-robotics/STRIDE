/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:36:11 GMT-05:00
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
  double t2469;
  double t2502;
  double t2539;
  double t2588;
  double t991;
  double t2650;
  double t2710;
  double t2714;
  double t2719;
  double t2725;
  double t2734;
  double t2741;
  double t2745;
  double t2679;
  double t2703;
  double t2709;
  double t2551;
  double t2610;
  double t2614;
  double t2659;
  double t2669;
  double t2672;
  double t2724;
  double t2727;
  double t2728;
  double t2753;
  double t2758;
  double t2759;
  double t2761;
  double t2810;
  double t2811;
  double t2823;
  t2469 = Cos(var1[5]);
  t2502 = -1.*t2469;
  t2539 = 1. + t2502;
  t2588 = Sin(var1[5]);
  t991 = Sin(var1[2]);
  t2650 = Cos(var1[2]);
  t2710 = Cos(var1[6]);
  t2714 = -1.*t2710;
  t2719 = 1. + t2714;
  t2725 = Sin(var1[6]);
  t2734 = t2650*t2469;
  t2741 = -1.*t991*t2588;
  t2745 = t2734 + t2741;
  t2679 = t2469*t991;
  t2703 = t2650*t2588;
  t2709 = t2679 + t2703;
  t2551 = -0.0695*t2539;
  t2610 = -0.0265*t2588;
  t2614 = t2551 + t2610;
  t2659 = -0.0265*t2539;
  t2669 = 0.0695*t2588;
  t2672 = t2659 + t2669;
  t2724 = -0.2375*t2719;
  t2727 = -0.0265*t2725;
  t2728 = t2724 + t2727;
  t2753 = -0.0265*t2719;
  t2758 = 0.2375*t2725;
  t2759 = t2753 + t2758;
  t2761 = t2710*t2745;
  t2810 = -1.*t2469*t991;
  t2811 = -1.*t2650*t2588;
  t2823 = t2810 + t2811;
  p_output1[0]=t2650*t2672 + t2709*t2728 - 0.2375*(t2709*t2710 + t2725*t2745) + t2745*t2759 - 0.0265*(-1.*t2709*t2725 + t2761) + t2614*t991 + var1[0];
  p_output1[1]=0.0293;
  p_output1[2]=t2614*t2650 + t2728*t2745 + t2759*t2823 - 0.0265*(-1.*t2725*t2745 + t2710*t2823) - 0.2375*(t2761 + t2725*t2823) - 1.*t2672*t991 + var1[1];
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
