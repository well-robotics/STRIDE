/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:07:20 GMT-05:00
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
  double t21;
  double t15;
  double t18;
  double t22;
  double t37;
  double t40;
  double t41;
  double t20;
  double t26;
  double t27;
  double t56;
  double t57;
  double t66;
  double t67;
  double t68;
  double t69;
  double t60;
  double t61;
  double t62;
  double t64;
  double t46;
  double t48;
  double t49;
  double t51;
  double t52;
  double t53;
  double t42;
  double t44;
  double t65;
  double t72;
  double t73;
  double t78;
  double t79;
  double t80;
  double t75;
  double t81;
  double t82;
  double t90;
  double t91;
  double t92;
  double t83;
  double t84;
  double t85;
  double t93;
  double t94;
  double t95;
  double t101;
  double t102;
  double t103;
  t21 = Cos(var1[2]);
  t15 = Cos(var1[5]);
  t18 = Sin(var1[2]);
  t22 = Sin(var1[5]);
  t37 = t21*t15;
  t40 = -1.*t18*t22;
  t41 = t37 + t40;
  t20 = t15*t18;
  t26 = t21*t22;
  t27 = t20 + t26;
  t56 = -1.*t15;
  t57 = 1. + t56;
  t66 = -0.0265*t57;
  t67 = -0.025413*t15;
  t68 = -0.08282*t22;
  t69 = t66 + t67 + t68;
  t60 = -0.0695*t57;
  t61 = -0.15232*t15;
  t62 = -0.0010869999999999977*t22;
  t64 = t60 + t61 + t62;
  t46 = -1.*t15*t18;
  t48 = -1.*t21*t22;
  t49 = t46 + t48;
  t51 = 0.69051*t49*t41;
  t52 = 0.69051*t27*t41;
  t53 = t51 + t52;
  t42 = Power(t41,2);
  t44 = 0.69051*t42;
  t65 = t15*t64;
  t72 = t69*t22;
  t73 = t65 + t72;
  t78 = -1.*t15*t69;
  t79 = t64*t22;
  t80 = t78 + t79;
  t75 = 0.69051*t41*t73;
  t81 = 0.69051*t27*t80;
  t82 = t75 + t81;
  t90 = 0.69051*t49*t73;
  t91 = 0.69051*t41*t80;
  t92 = t90 + t91;
  t83 = -0.0007505843699999984*t27;
  t84 = -0.0571880382*t41;
  t85 = t83 + t84;
  t93 = -0.0571880382*t49;
  t94 = -0.0007505843699999984*t41;
  t95 = t93 + t94;
  t101 = -0.0571880382*t73;
  t102 = -0.0007505843699999984*t80;
  t103 = 0.00025 + t101 + t102;
  p_output1[0]=0.69051*Power(t27,2) + t44;
  p_output1[1]=t53;
  p_output1[2]=t82;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t85;
  p_output1[6]=0;
  p_output1[7]=t53;
  p_output1[8]=t44 + 0.69051*Power(t49,2);
  p_output1[9]=t92;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t95;
  p_output1[13]=0;
  p_output1[14]=t82;
  p_output1[15]=t92;
  p_output1[16]=0.00025 + 0.69051*Power(t73,2) + 0.69051*Power(t80,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t103;
  p_output1[20]=0;
  p_output1[21]=0;
  p_output1[22]=0;
  p_output1[23]=0;
  p_output1[24]=0;
  p_output1[25]=0;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=0;
  p_output1[29]=0;
  p_output1[30]=0;
  p_output1[31]=0;
  p_output1[32]=0;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=t85;
  p_output1[36]=t95;
  p_output1[37]=t103;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.004987129208934191;
  p_output1[41]=0;
  p_output1[42]=0;
  p_output1[43]=0;
  p_output1[44]=0;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=0;
  p_output1[48]=0;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "Mmat4_Robot_Assembly_v3_straight_leg_TPU_foot.hh"

namespace SymFunction
{

void Mmat4_Robot_Assembly_v3_straight_leg_TPU_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
