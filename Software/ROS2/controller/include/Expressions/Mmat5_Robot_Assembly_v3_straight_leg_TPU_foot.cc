/*
 * Automatically Generated from Mathematica.
 * Sun 9 Jun 2024 16:07:21 GMT-05:00
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
  double t31;
  double t26;
  double t27;
  double t37;
  double t45;
  double t20;
  double t49;
  double t51;
  double t52;
  double t30;
  double t41;
  double t42;
  double t44;
  double t53;
  double t60;
  double t64;
  double t65;
  double t66;
  double t67;
  double t68;
  double t72;
  double t93;
  double t94;
  double t89;
  double t90;
  double t95;
  double t96;
  double t97;
  double t98;
  double t100;
  double t101;
  double t102;
  double t103;
  double t107;
  double t108;
  double t109;
  double t110;
  double t111;
  double t91;
  double t92;
  double t99;
  double t104;
  double t105;
  double t80;
  double t81;
  double t82;
  double t83;
  double t84;
  double t85;
  double t87;
  double t88;
  double t106;
  double t112;
  double t113;
  double t115;
  double t116;
  double t117;
  double t126;
  double t127;
  double t128;
  double t129;
  double t130;
  double t120;
  double t121;
  double t122;
  double t123;
  double t124;
  double t114;
  double t118;
  double t119;
  double t141;
  double t142;
  double t143;
  double t125;
  double t131;
  double t132;
  double t144;
  double t145;
  double t146;
  double t155;
  double t156;
  double t157;
  double t133;
  double t134;
  double t135;
  double t147;
  double t148;
  double t149;
  double t158;
  double t159;
  double t160;
  double t166;
  double t167;
  double t168;
  t31 = Cos(var1[5]);
  t26 = Cos(var1[6]);
  t27 = Sin(var1[5]);
  t37 = Sin(var1[6]);
  t45 = Cos(var1[2]);
  t20 = Sin(var1[2]);
  t49 = t31*t26;
  t51 = -1.*t27*t37;
  t52 = t49 + t51;
  t30 = -1.*t26*t27;
  t41 = -1.*t31*t37;
  t42 = t30 + t41;
  t44 = t20*t42;
  t53 = t45*t52;
  t60 = t44 + t53;
  t64 = t26*t27;
  t65 = t31*t37;
  t66 = t64 + t65;
  t67 = t45*t66;
  t68 = t20*t52;
  t72 = t67 + t68;
  t93 = -1.*t26;
  t94 = 1. + t93;
  t89 = -1.*t31;
  t90 = 1. + t89;
  t95 = -0.0265*t94;
  t96 = -0.025157*t26;
  t97 = -0.07483400000000001*t37;
  t98 = t95 + t96 + t97;
  t100 = -0.2375*t94;
  t101 = -0.312334*t26;
  t102 = -0.0013430000000000004*t37;
  t103 = t100 + t101 + t102;
  t107 = -0.0695*t90;
  t108 = -0.0265*t27;
  t109 = -1.*t27*t98;
  t110 = t31*t103;
  t111 = t107 + t108 + t109 + t110;
  t91 = -0.0265*t90;
  t92 = 0.0695*t27;
  t99 = t31*t98;
  t104 = t27*t103;
  t105 = t91 + t92 + t99 + t104;
  t80 = t45*t42;
  t81 = -1.*t20*t52;
  t82 = t80 + t81;
  t83 = 0.19605*t60*t82;
  t84 = -1.*t20*t66;
  t85 = t84 + t53;
  t87 = 0.19605*t85*t72;
  t88 = t83 + t87;
  t106 = t105*t66;
  t112 = t111*t52;
  t113 = t106 + t112;
  t115 = -1.*t111*t42;
  t116 = -1.*t105*t52;
  t117 = t115 + t116;
  t126 = -0.0265*t26;
  t127 = -1.*t26*t98;
  t128 = 0.0695*t37;
  t129 = t103*t37;
  t130 = t126 + t127 + t128 + t129;
  t120 = 0.0695*t26;
  t121 = t26*t103;
  t122 = 0.0265*t37;
  t123 = t98*t37;
  t124 = t120 + t121 + t122 + t123;
  t114 = 0.19605*t60*t113;
  t118 = 0.19605*t72*t117;
  t119 = t114 + t118;
  t141 = 0.19605*t82*t113;
  t142 = 0.19605*t85*t117;
  t143 = t141 + t142;
  t125 = 0.19605*t124*t60;
  t131 = 0.19605*t130*t72;
  t132 = t125 + t131;
  t144 = 0.19605*t130*t85;
  t145 = 0.19605*t124*t82;
  t146 = t144 + t145;
  t155 = 0.19605*t124*t113;
  t156 = 0.19605*t130*t117;
  t157 = 0.00011 + t155 + t156;
  t133 = -0.014671205700000002*t60;
  t134 = -0.0002632951500000001*t72;
  t135 = t133 + t134;
  t147 = -0.0002632951500000001*t85;
  t148 = -0.014671205700000002*t82;
  t149 = t147 + t148;
  t158 = -0.014671205700000002*t113;
  t159 = -0.0002632951500000001*t117;
  t160 = 0.00011 + t158 + t159;
  t166 = -0.014671205700000002*t124;
  t167 = -0.0002632951500000001*t130;
  t168 = 0.00011 + t166 + t167;
  p_output1[0]=0.19605*Power(t60,2) + 0.19605*Power(t72,2);
  p_output1[1]=t88;
  p_output1[2]=t119;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t132;
  p_output1[6]=t135;
  p_output1[7]=t88;
  p_output1[8]=0.19605*Power(t82,2) + 0.19605*Power(t85,2);
  p_output1[9]=t143;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t146;
  p_output1[13]=t149;
  p_output1[14]=t119;
  p_output1[15]=t143;
  p_output1[16]=0.00011 + 0.19605*Power(t113,2) + 0.19605*Power(t117,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t157;
  p_output1[20]=t160;
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
  p_output1[35]=t132;
  p_output1[36]=t146;
  p_output1[37]=t157;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.00011 + 0.19605*Power(t124,2) + 0.19605*Power(t130,2);
  p_output1[41]=t168;
  p_output1[42]=t135;
  p_output1[43]=t149;
  p_output1[44]=t160;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=t168;
  p_output1[48]=0.0012082586127402503;
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

#include "Mmat5_Robot_Assembly_v3_straight_leg_TPU_foot.hh"

namespace SymFunction
{

void Mmat5_Robot_Assembly_v3_straight_leg_TPU_foot_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
