/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:27 GMT-05:00
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
  double t66;
  double t73;
  double t76;
  double t82;
  double t62;
  double t123;
  double t155;
  double t157;
  double t158;
  double t164;
  double t148;
  double t150;
  double t153;
  double t170;
  double t171;
  double t172;
  double t78;
  double t92;
  double t101;
  double t125;
  double t133;
  double t136;
  double t161;
  double t166;
  double t167;
  double t175;
  double t186;
  double t191;
  double t201;
  double t216;
  double t217;
  double t219;
  t66 = Cos(var1[3]);
  t73 = -1.*t66;
  t76 = 1. + t73;
  t82 = Sin(var1[3]);
  t62 = Cos(var1[2]);
  t123 = Sin(var1[2]);
  t155 = Cos(var1[4]);
  t157 = -1.*t155;
  t158 = 1. + t157;
  t164 = Sin(var1[4]);
  t148 = t62*t66;
  t150 = t123*t82;
  t153 = t148 + t150;
  t170 = t66*t123;
  t171 = -1.*t62*t82;
  t172 = t170 + t171;
  t78 = -0.0265*t76;
  t92 = -0.0695*t82;
  t101 = t78 + t92;
  t125 = -0.0695*t76;
  t133 = 0.0265*t82;
  t136 = t125 + t133;
  t161 = -0.0265*t158;
  t166 = -0.2375*t164;
  t167 = t161 + t166;
  t175 = -0.2375*t158;
  t186 = 0.0265*t164;
  t191 = t175 + t186;
  t201 = t155*t153;
  t216 = -1.*t66*t123;
  t217 = t62*t82;
  t219 = t216 + t217;
  p_output1[0]=t123*t136 + t153*t167 + 0.0225*(-1.*t153*t164 + t155*t172) + t172*t191 - 0.0265*(t164*t172 + t201) + t101*t62 + var1[0];
  p_output1[1]=-1.*t101*t123 + t153*t191 + t167*t219 - 0.0265*(t153*t164 + t155*t219) + 0.0225*(t201 - 1.*t164*t219) + t136*t62 + var1[1];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "pRightToe.hh"

namespace SymFunction
{

void pRightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
