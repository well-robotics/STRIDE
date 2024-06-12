/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:36:46 GMT-05:00
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
  double t49;
  double t39;
  double t43;
  double t51;
  double t65;
  double t38;
  double t66;
  double t67;
  double t70;
  double t48;
  double t53;
  double t59;
  double t60;
  double t71;
  double t72;
  double t78;
  double t80;
  double t81;
  double t83;
  double t84;
  double t86;
  double t111;
  double t112;
  double t114;
  double t116;
  double t117;
  double t118;
  double t123;
  double t124;
  double t125;
  double t131;
  double t133;
  double t134;
  double t115;
  double t126;
  double t127;
  double t96;
  double t98;
  double t99;
  double t101;
  double t103;
  double t105;
  double t106;
  double t107;
  double t129;
  double t135;
  double t138;
  double t140;
  double t141;
  double t142;
  double t151;
  double t152;
  double t153;
  double t145;
  double t147;
  double t148;
  double t139;
  double t143;
  double t144;
  double t170;
  double t171;
  double t172;
  double t149;
  double t154;
  double t155;
  double t175;
  double t176;
  double t177;
  double t187;
  double t188;
  double t189;
  double t157;
  double t160;
  double t161;
  double t178;
  double t179;
  double t180;
  double t190;
  double t191;
  double t192;
  double t199;
  double t200;
  double t201;
  t49 = Cos(var1[5]);
  t39 = Cos(var1[6]);
  t43 = Sin(var1[5]);
  t51 = Sin(var1[6]);
  t65 = Cos(var1[2]);
  t38 = Sin(var1[2]);
  t66 = t49*t39;
  t67 = -1.*t43*t51;
  t70 = t66 + t67;
  t48 = -1.*t39*t43;
  t53 = -1.*t49*t51;
  t59 = t48 + t53;
  t60 = t38*t59;
  t71 = t65*t70;
  t72 = t60 + t71;
  t78 = t39*t43;
  t80 = t49*t51;
  t81 = t78 + t80;
  t83 = t65*t81;
  t84 = t38*t70;
  t86 = t83 + t84;
  t111 = -0.022659*t39;
  t112 = -0.007367999999999986*t51;
  t114 = t111 + t112;
  t116 = -1.*t39;
  t117 = 1. + t116;
  t118 = -0.16*t117;
  t123 = -0.167368*t39;
  t124 = 0.022659*t51;
  t125 = t118 + t123 + t124;
  t131 = -1.*t43*t114;
  t133 = t49*t125;
  t134 = t131 + t133;
  t115 = t49*t114;
  t126 = t43*t125;
  t127 = t115 + t126;
  t96 = t65*t59;
  t98 = -1.*t38*t70;
  t99 = t96 + t98;
  t101 = 0.14994*t72*t99;
  t103 = -1.*t38*t81;
  t105 = t103 + t71;
  t106 = 0.14994*t105*t86;
  t107 = t101 + t106;
  t129 = t127*t81;
  t135 = t134*t70;
  t138 = t129 + t135;
  t140 = -1.*t134*t59;
  t141 = -1.*t127*t70;
  t142 = t140 + t141;
  t151 = -1.*t39*t114;
  t152 = t125*t51;
  t153 = t151 + t152;
  t145 = t39*t125;
  t147 = t114*t51;
  t148 = t145 + t147;
  t139 = 0.14994*t72*t138;
  t143 = 0.14994*t86*t142;
  t144 = t139 + t143;
  t170 = 0.14994*t99*t138;
  t171 = 0.14994*t105*t142;
  t172 = t170 + t171;
  t149 = 0.14994*t148*t72;
  t154 = 0.14994*t153*t86;
  t155 = t149 + t154;
  t175 = 0.14994*t153*t105;
  t176 = 0.14994*t148*t99;
  t177 = t175 + t176;
  t187 = 0.14994*t148*t138;
  t188 = 0.14994*t153*t142;
  t189 = 0.000088 + t187 + t188;
  t157 = -0.0011047579199999977*t72;
  t160 = 0.0033974904599999994*t86;
  t161 = t157 + t160;
  t178 = 0.0033974904599999994*t105;
  t179 = -0.0011047579199999977*t99;
  t180 = t178 + t179;
  t190 = -0.0011047579199999977*t138;
  t191 = 0.0033974904599999994*t142;
  t192 = 0.000088 + t190 + t191;
  t199 = -0.0011047579199999977*t148;
  t200 = 0.0033974904599999994*t153;
  t201 = 0.000088 + t199 + t200;
  p_output1[0]=0.14994*Power(t72,2) + 0.14994*Power(t86,2);
  p_output1[1]=t107;
  p_output1[2]=t144;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t155;
  p_output1[6]=t161;
  p_output1[7]=t107;
  p_output1[8]=0.14994*Power(t105,2) + 0.14994*Power(t99,2);
  p_output1[9]=t172;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t177;
  p_output1[13]=t180;
  p_output1[14]=t144;
  p_output1[15]=t172;
  p_output1[16]=0.000088 + 0.14994*Power(t138,2) + 0.14994*Power(t142,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t189;
  p_output1[20]=t192;
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
  p_output1[35]=t155;
  p_output1[36]=t177;
  p_output1[37]=t189;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.000088 + 0.14994*Power(t148,2) + 0.14994*Power(t153,2);
  p_output1[41]=t201;
  p_output1[42]=t161;
  p_output1[43]=t180;
  p_output1[44]=t192;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=t201;
  p_output1[48]=0.00017312359268769993;
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

#include "Mmat5_five_link_walker.hh"

namespace SymFunction
{

void Mmat5_five_link_walker_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
