/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:21 GMT-05:00
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
  double t5;
  double t7;
  double t13;
  double t14;
  double t15;
  double t17;
  double t16;
  double t18;
  double t19;
  double t20;
  double t21;
  double t22;
  double t23;
  double t24;
  double t29;
  double t30;
  double t31;
  double t25;
  double t26;
  double t27;
  double t35;
  double t36;
  double t37;
  double t39;
  double t57;
  double t58;
  double t59;
  double t61;
  double t60;
  double t62;
  double t63;
  double t64;
  double t65;
  double t66;
  double t67;
  double t68;
  double t69;
  double t70;
  double t71;
  double t73;
  double t74;
  double t75;
  double t79;
  double t80;
  double t81;
  double t83;
  double t107;
  double t108;
  double t109;
  double t110;
  double t111;
  double t38;
  double t40;
  double t41;
  double t43;
  double t44;
  double t45;
  double t47;
  double t127;
  double t128;
  double t82;
  double t84;
  double t85;
  double t129;
  double t130;
  double t131;
  double t87;
  double t88;
  double t89;
  double t91;
  t5 = Cos(var1[2]);
  t7 = Sin(var1[2]);
  t13 = Cos(var1[3]);
  t14 = -1.*t13;
  t15 = 1. + t14;
  t17 = Sin(var1[3]);
  t16 = -0.0265*t15;
  t18 = -0.0695*t17;
  t19 = t16 + t18;
  t20 = t5*t19;
  t21 = -0.0695*t15;
  t22 = 0.0265*t17;
  t23 = t21 + t22;
  t24 = t7*t23;
  t29 = t5*t13;
  t30 = t7*t17;
  t31 = t29 + t30;
  t25 = t13*t7;
  t26 = -1.*t5*t17;
  t27 = t25 + t26;
  t35 = Cos(var1[4]);
  t36 = -1.*t35;
  t37 = 1. + t36;
  t39 = Sin(var1[4]);
  t57 = Cos(var1[5]);
  t58 = -1.*t57;
  t59 = 1. + t58;
  t61 = Sin(var1[5]);
  t60 = -0.0695*t59;
  t62 = -0.0265*t61;
  t63 = t60 + t62;
  t64 = t7*t63;
  t65 = -0.0265*t59;
  t66 = 0.0695*t61;
  t67 = t65 + t66;
  t68 = t5*t67;
  t69 = t57*t7;
  t70 = t5*t61;
  t71 = t69 + t70;
  t73 = t5*t57;
  t74 = -1.*t7*t61;
  t75 = t73 + t74;
  t79 = Cos(var1[6]);
  t80 = -1.*t79;
  t81 = 1. + t80;
  t83 = Sin(var1[6]);
  t107 = -1.*t7*t19;
  t108 = t5*t23;
  t109 = -1.*t13*t7;
  t110 = t5*t17;
  t111 = t109 + t110;
  t38 = -0.0265*t37;
  t40 = -0.2375*t39;
  t41 = t38 + t40;
  t43 = -0.2375*t37;
  t44 = 0.0265*t39;
  t45 = t43 + t44;
  t47 = t35*t31;
  t127 = t5*t63;
  t128 = -1.*t7*t67;
  t82 = -0.2375*t81;
  t84 = -0.0265*t83;
  t85 = t82 + t84;
  t129 = -1.*t57*t7;
  t130 = -1.*t5*t61;
  t131 = t129 + t130;
  t87 = -0.0265*t81;
  t88 = 0.2375*t83;
  t89 = t87 + t88;
  t91 = t79*t75;
  p_output1[0]=0.2996793431028799*(0.69051*(t20 + t24 - 0.15232*t27 - 0.025367*t31 + var1[0]) + 0.19964*(t20 + t24 - 0.314514*(t27*t35 - 1.*t31*t39) + t31*t41 + t27*t45 - 0.025229*(t27*t39 + t47) + var1[0]) + 1.5566*(-0.026461*t5 - 0.046589*t7 + var1[0]) + 0.69051*(t64 + t68 - 0.15232*t71 - 0.025413*t75 + var1[0]) + 0.19964*(t64 + t68 - 0.314506*(t71*t79 + t75*t83) + t71*t85 + t75*t89 - 0.025226*(-1.*t71*t83 + t91) + var1[0]));
  p_output1[1]=0.002039417585183853;
  p_output1[2]=0.2996793431028799*(0.69051*(t107 + t108 - 0.025367*t111 - 0.15232*t31 + var1[1]) + 0.19964*(t107 + t108 - 0.025229*(t111*t35 + t31*t39) + t111*t41 + t31*t45 - 0.314514*(-1.*t111*t39 + t47) + var1[1]) + 1.5566*(-0.046589*t5 + 0.026461*t7 + var1[1]) + 0.69051*(t127 + t128 - 0.025413*t131 - 0.15232*t75 + var1[1]) + 0.19964*(t127 + t128 - 0.025226*(t131*t79 - 1.*t75*t83) + t75*t85 + t131*t89 - 0.314506*(t131*t83 + t91) + var1[1]));
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

#include "COMPosition.hh"

namespace SymFunction
{

void COMPosition_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
