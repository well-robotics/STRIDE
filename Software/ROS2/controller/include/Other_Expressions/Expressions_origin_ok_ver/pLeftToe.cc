/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:26 GMT-05:00
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
  double t46;
  double t55;
  double t60;
  double t64;
  double t42;
  double t77;
  double t116;
  double t119;
  double t120;
  double t124;
  double t137;
  double t138;
  double t139;
  double t101;
  double t113;
  double t115;
  double t62;
  double t66;
  double t73;
  double t78;
  double t82;
  double t89;
  double t123;
  double t125;
  double t133;
  double t148;
  double t150;
  double t153;
  double t156;
  double t176;
  double t178;
  double t183;
  t46 = Cos(var1[5]);
  t55 = -1.*t46;
  t60 = 1. + t55;
  t64 = Sin(var1[5]);
  t42 = Sin(var1[2]);
  t77 = Cos(var1[2]);
  t116 = Cos(var1[6]);
  t119 = -1.*t116;
  t120 = 1. + t119;
  t124 = Sin(var1[6]);
  t137 = t77*t46;
  t138 = -1.*t42*t64;
  t139 = t137 + t138;
  t101 = t46*t42;
  t113 = t77*t64;
  t115 = t101 + t113;
  t62 = -0.0695*t60;
  t66 = -0.0265*t64;
  t73 = t62 + t66;
  t78 = -0.0265*t60;
  t82 = 0.0695*t64;
  t89 = t78 + t82;
  t123 = -0.2375*t120;
  t125 = -0.0265*t124;
  t133 = t123 + t125;
  t148 = -0.0265*t120;
  t150 = 0.2375*t124;
  t153 = t148 + t150;
  t156 = t116*t139;
  t176 = -1.*t46*t42;
  t178 = -1.*t77*t64;
  t183 = t176 + t178;
  p_output1[0]=t115*t133 + 0.0225*(t115*t116 + t124*t139) + t139*t153 - 0.0265*(-1.*t115*t124 + t156) + t42*t73 + t77*t89 + var1[0];
  p_output1[1]=t133*t139 + t153*t183 - 0.0265*(-1.*t124*t139 + t116*t183) + 0.0225*(t156 + t124*t183) + t73*t77 - 1.*t42*t89 + var1[1];
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

#include "pLeftToe.hh"

namespace SymFunction
{

void pLeftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
