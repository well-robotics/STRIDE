/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:10 GMT-05:00
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
  double t26;
  double t17;
  double t20;
  double t27;
  double t49;
  double t53;
  double t59;
  double t24;
  double t33;
  double t37;
  double t78;
  double t79;
  double t92;
  double t93;
  double t95;
  double t96;
  double t83;
  double t84;
  double t86;
  double t90;
  double t66;
  double t68;
  double t70;
  double t72;
  double t73;
  double t74;
  double t60;
  double t64;
  double t91;
  double t100;
  double t101;
  double t110;
  double t111;
  double t112;
  double t105;
  double t115;
  double t116;
  double t128;
  double t129;
  double t130;
  double t117;
  double t119;
  double t120;
  double t131;
  double t132;
  double t133;
  double t144;
  double t145;
  double t146;
  t26 = Cos(var1[2]);
  t17 = Cos(var1[5]);
  t20 = Sin(var1[2]);
  t27 = Sin(var1[5]);
  t49 = t26*t17;
  t53 = -1.*t20*t27;
  t59 = t49 + t53;
  t24 = t17*t20;
  t33 = t26*t27;
  t37 = t24 + t33;
  t78 = -1.*t17;
  t79 = 1. + t78;
  t92 = -0.0265*t79;
  t93 = -0.025413*t17;
  t95 = -0.08282*t27;
  t96 = t92 + t93 + t95;
  t83 = -0.0695*t79;
  t84 = -0.15232*t17;
  t86 = -0.0010869999999999977*t27;
  t90 = t83 + t84 + t86;
  t66 = -1.*t17*t20;
  t68 = -1.*t26*t27;
  t70 = t66 + t68;
  t72 = 0.69051*t70*t59;
  t73 = 0.69051*t37*t59;
  t74 = t72 + t73;
  t60 = Power(t59,2);
  t64 = 0.69051*t60;
  t91 = t17*t90;
  t100 = t96*t27;
  t101 = t91 + t100;
  t110 = -1.*t17*t96;
  t111 = t90*t27;
  t112 = t110 + t111;
  t105 = 0.69051*t59*t101;
  t115 = 0.69051*t37*t112;
  t116 = t105 + t115;
  t128 = 0.69051*t70*t101;
  t129 = 0.69051*t59*t112;
  t130 = t128 + t129;
  t117 = -0.0007505843699999984*t37;
  t119 = -0.0571880382*t59;
  t120 = t117 + t119;
  t131 = -0.0571880382*t70;
  t132 = -0.0007505843699999984*t59;
  t133 = t131 + t132;
  t144 = -0.0571880382*t101;
  t145 = -0.0007505843699999984*t112;
  t146 = 0.00025 + t144 + t145;
  p_output1[0]=0.69051*Power(t37,2) + t64;
  p_output1[1]=t74;
  p_output1[2]=t116;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t120;
  p_output1[6]=0;
  p_output1[7]=t74;
  p_output1[8]=t64 + 0.69051*Power(t70,2);
  p_output1[9]=t130;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t133;
  p_output1[13]=0;
  p_output1[14]=t116;
  p_output1[15]=t130;
  p_output1[16]=0.00025 + 0.69051*Power(t101,2) + 0.69051*Power(t112,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t146;
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
  p_output1[35]=t120;
  p_output1[36]=t133;
  p_output1[37]=t146;
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

#include "Mmat4_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Mmat4_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
