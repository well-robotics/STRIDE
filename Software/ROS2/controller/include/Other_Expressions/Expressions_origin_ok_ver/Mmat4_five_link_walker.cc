/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:36:44 GMT-05:00
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
  double t30;
  double t19;
  double t26;
  double t31;
  double t39;
  double t15;
  double t43;
  double t45;
  double t46;
  double t27;
  double t35;
  double t37;
  double t38;
  double t48;
  double t49;
  double t59;
  double t60;
  double t63;
  double t65;
  double t66;
  double t67;
  double t84;
  double t86;
  double t87;
  double t91;
  double t92;
  double t93;
  double t95;
  double t96;
  double t98;
  double t102;
  double t103;
  double t104;
  double t90;
  double t99;
  double t100;
  double t73;
  double t74;
  double t75;
  double t78;
  double t79;
  double t80;
  double t81;
  double t83;
  double t101;
  double t105;
  double t106;
  double t108;
  double t110;
  double t111;
  double t120;
  double t123;
  double t124;
  double t115;
  double t116;
  double t117;
  double t107;
  double t112;
  double t114;
  double t140;
  double t141;
  double t142;
  double t118;
  double t125;
  double t126;
  double t143;
  double t144;
  double t145;
  double t157;
  double t160;
  double t161;
  double t127;
  double t129;
  double t131;
  double t147;
  double t148;
  double t149;
  double t162;
  double t163;
  double t164;
  double t175;
  double t176;
  double t177;
  t30 = Cos(var1[3]);
  t19 = Cos(var1[4]);
  t26 = Sin(var1[3]);
  t31 = Sin(var1[4]);
  t39 = Cos(var1[2]);
  t15 = Sin(var1[2]);
  t43 = t30*t19;
  t45 = -1.*t26*t31;
  t46 = t43 + t45;
  t27 = -1.*t19*t26;
  t35 = -1.*t30*t31;
  t37 = t27 + t35;
  t38 = t15*t37;
  t48 = t39*t46;
  t49 = t38 + t48;
  t59 = t19*t26;
  t60 = t30*t31;
  t63 = t59 + t60;
  t65 = t39*t63;
  t66 = t15*t46;
  t67 = t65 + t66;
  t84 = -0.022663*t19;
  t86 = -0.007370999999999989*t31;
  t87 = t84 + t86;
  t91 = -1.*t19;
  t92 = 1. + t91;
  t93 = -0.16*t92;
  t95 = -0.167371*t19;
  t96 = 0.022663*t31;
  t98 = t93 + t95 + t96;
  t102 = -1.*t26*t87;
  t103 = t30*t98;
  t104 = t102 + t103;
  t90 = t30*t87;
  t99 = t26*t98;
  t100 = t90 + t99;
  t73 = t39*t37;
  t74 = -1.*t15*t46;
  t75 = t73 + t74;
  t78 = 0.14994*t49*t75;
  t79 = -1.*t15*t63;
  t80 = t79 + t48;
  t81 = 0.14994*t80*t67;
  t83 = t78 + t81;
  t101 = t100*t63;
  t105 = t104*t46;
  t106 = t101 + t105;
  t108 = -1.*t104*t37;
  t110 = -1.*t100*t46;
  t111 = t108 + t110;
  t120 = -1.*t19*t87;
  t123 = t98*t31;
  t124 = t120 + t123;
  t115 = t19*t98;
  t116 = t87*t31;
  t117 = t115 + t116;
  t107 = 0.14994*t49*t106;
  t112 = 0.14994*t67*t111;
  t114 = t107 + t112;
  t140 = 0.14994*t75*t106;
  t141 = 0.14994*t80*t111;
  t142 = t140 + t141;
  t118 = 0.14994*t117*t49;
  t125 = 0.14994*t124*t67;
  t126 = t118 + t125;
  t143 = 0.14994*t124*t80;
  t144 = 0.14994*t117*t75;
  t145 = t143 + t144;
  t157 = 0.14994*t117*t106;
  t160 = 0.14994*t124*t111;
  t161 = 0.000088 + t157 + t160;
  t127 = -0.0011052077399999983*t49;
  t129 = 0.0033980902199999994*t67;
  t131 = t127 + t129;
  t147 = 0.0033980902199999994*t80;
  t148 = -0.0011052077399999983*t75;
  t149 = t147 + t148;
  t162 = -0.0011052077399999983*t106;
  t163 = 0.0033980902199999994*t111;
  t164 = 0.000088 + t162 + t163;
  t175 = -0.0011052077399999983*t117;
  t176 = 0.0033980902199999994*t124;
  t177 = 0.000088 + t175 + t176;
  p_output1[0]=0.14994*Power(t49,2) + 0.14994*Power(t67,2);
  p_output1[1]=t83;
  p_output1[2]=t114;
  p_output1[3]=t126;
  p_output1[4]=t131;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t83;
  p_output1[8]=0.14994*Power(t75,2) + 0.14994*Power(t80,2);
  p_output1[9]=t142;
  p_output1[10]=t145;
  p_output1[11]=t149;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t114;
  p_output1[15]=t142;
  p_output1[16]=0.000088 + 0.14994*Power(t106,2) + 0.14994*Power(t111,2);
  p_output1[17]=t161;
  p_output1[18]=t164;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t126;
  p_output1[22]=t145;
  p_output1[23]=t161;
  p_output1[24]=0.000088 + 0.14994*Power(t117,2) + 0.14994*Power(t124,2);
  p_output1[25]=t177;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=t131;
  p_output1[29]=t149;
  p_output1[30]=t164;
  p_output1[31]=t177;
  p_output1[32]=0.00017315740490739996;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=0;
  p_output1[36]=0;
  p_output1[37]=0;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0;
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

#include "Mmat4_five_link_walker.hh"

namespace SymFunction
{

void Mmat4_five_link_walker_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
