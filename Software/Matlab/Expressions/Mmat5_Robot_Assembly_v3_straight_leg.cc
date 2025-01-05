/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:13 GMT-05:00
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
  double t41;
  double t33;
  double t37;
  double t49;
  double t65;
  double t24;
  double t70;
  double t72;
  double t73;
  double t40;
  double t59;
  double t60;
  double t64;
  double t74;
  double t83;
  double t90;
  double t91;
  double t92;
  double t93;
  double t95;
  double t100;
  double t131;
  double t132;
  double t126;
  double t128;
  double t133;
  double t136;
  double t138;
  double t141;
  double t143;
  double t144;
  double t145;
  double t146;
  double t152;
  double t153;
  double t154;
  double t155;
  double t156;
  double t129;
  double t130;
  double t142;
  double t149;
  double t150;
  double t112;
  double t115;
  double t116;
  double t117;
  double t119;
  double t120;
  double t124;
  double t125;
  double t151;
  double t158;
  double t159;
  double t162;
  double t163;
  double t164;
  double t176;
  double t177;
  double t179;
  double t180;
  double t182;
  double t169;
  double t170;
  double t172;
  double t173;
  double t174;
  double t160;
  double t165;
  double t168;
  double t197;
  double t198;
  double t199;
  double t175;
  double t184;
  double t185;
  double t200;
  double t201;
  double t202;
  double t211;
  double t212;
  double t213;
  double t186;
  double t189;
  double t190;
  double t203;
  double t204;
  double t205;
  double t214;
  double t215;
  double t216;
  double t223;
  double t224;
  double t225;
  t41 = Cos(var1[5]);
  t33 = Cos(var1[6]);
  t37 = Sin(var1[5]);
  t49 = Sin(var1[6]);
  t65 = Cos(var1[2]);
  t24 = Sin(var1[2]);
  t70 = t41*t33;
  t72 = -1.*t37*t49;
  t73 = t70 + t72;
  t40 = -1.*t33*t37;
  t59 = -1.*t41*t49;
  t60 = t40 + t59;
  t64 = t24*t60;
  t74 = t65*t73;
  t83 = t64 + t74;
  t90 = t33*t37;
  t91 = t41*t49;
  t92 = t90 + t91;
  t93 = t65*t92;
  t95 = t24*t73;
  t100 = t93 + t95;
  t131 = -1.*t33;
  t132 = 1. + t131;
  t126 = -1.*t41;
  t128 = 1. + t126;
  t133 = -0.0265*t132;
  t136 = -0.025226*t33;
  t138 = -0.07700600000000002*t49;
  t141 = t133 + t136 + t138;
  t143 = -0.2375*t132;
  t144 = -0.314506*t33;
  t145 = -0.0012740000000000008*t49;
  t146 = t143 + t144 + t145;
  t152 = -0.0695*t128;
  t153 = -0.0265*t37;
  t154 = -1.*t37*t141;
  t155 = t41*t146;
  t156 = t152 + t153 + t154 + t155;
  t129 = -0.0265*t128;
  t130 = 0.0695*t37;
  t142 = t41*t141;
  t149 = t37*t146;
  t150 = t129 + t130 + t142 + t149;
  t112 = t65*t60;
  t115 = -1.*t24*t73;
  t116 = t112 + t115;
  t117 = 0.19964*t83*t116;
  t119 = -1.*t24*t92;
  t120 = t119 + t74;
  t124 = 0.19964*t120*t100;
  t125 = t117 + t124;
  t151 = t150*t92;
  t158 = t156*t73;
  t159 = t151 + t158;
  t162 = -1.*t156*t60;
  t163 = -1.*t150*t73;
  t164 = t162 + t163;
  t176 = -0.0265*t33;
  t177 = -1.*t33*t141;
  t179 = 0.0695*t49;
  t180 = t146*t49;
  t182 = t176 + t177 + t179 + t180;
  t169 = 0.0695*t33;
  t170 = t33*t146;
  t172 = 0.0265*t49;
  t173 = t141*t49;
  t174 = t169 + t170 + t172 + t173;
  t160 = 0.19964*t83*t159;
  t165 = 0.19964*t100*t164;
  t168 = t160 + t165;
  t197 = 0.19964*t116*t159;
  t198 = 0.19964*t120*t164;
  t199 = t197 + t198;
  t175 = 0.19964*t174*t83;
  t184 = 0.19964*t182*t100;
  t185 = t175 + t184;
  t200 = 0.19964*t182*t120;
  t201 = 0.19964*t174*t116;
  t202 = t200 + t201;
  t211 = 0.19964*t174*t159;
  t212 = 0.19964*t182*t164;
  t213 = 0.000105 + t211 + t212;
  t186 = -0.015373477840000005*t83;
  t189 = -0.0002543413600000002*t100;
  t190 = t186 + t189;
  t203 = -0.0002543413600000002*t120;
  t204 = -0.015373477840000005*t116;
  t205 = t203 + t204;
  t214 = -0.015373477840000005*t159;
  t215 = -0.0002543413600000002*t164;
  t216 = 0.000105 + t214 + t215;
  t223 = -0.015373477840000005*t174;
  t224 = -0.0002543413600000002*t182;
  t225 = 0.000105 + t223 + t224;
  p_output1[0]=0.19964*Power(t100,2) + 0.19964*Power(t83,2);
  p_output1[1]=t125;
  p_output1[2]=t168;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t185;
  p_output1[6]=t190;
  p_output1[7]=t125;
  p_output1[8]=0.19964*Power(t116,2) + 0.19964*Power(t120,2);
  p_output1[9]=t199;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t202;
  p_output1[13]=t205;
  p_output1[14]=t168;
  p_output1[15]=t199;
  p_output1[16]=0.000105 + 0.19964*Power(t159,2) + 0.19964*Power(t164,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t213;
  p_output1[20]=t216;
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
  p_output1[35]=t185;
  p_output1[36]=t202;
  p_output1[37]=t213;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.000105 + 0.19964*Power(t174,2) + 0.19964*Power(t182,2);
  p_output1[41]=t225;
  p_output1[42]=t190;
  p_output1[43]=t205;
  p_output1[44]=t216;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=t225;
  p_output1[48]=0.0012891740654396807;
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

#include "Mmat5_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Mmat5_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
