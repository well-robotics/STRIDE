/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:22 GMT-05:00
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
  double t9;
  double t6;
  double t20;
  double t21;
  double t22;
  double t24;
  double t23;
  double t25;
  double t28;
  double t29;
  double t30;
  double t32;
  double t33;
  double t34;
  double t38;
  double t40;
  double t41;
  double t43;
  double t44;
  double t45;
  double t50;
  double t51;
  double t52;
  double t54;
  double t84;
  double t85;
  double t86;
  double t88;
  double t87;
  double t89;
  double t90;
  double t92;
  double t93;
  double t94;
  double t95;
  double t96;
  double t101;
  double t102;
  double t103;
  double t97;
  double t98;
  double t99;
  double t108;
  double t110;
  double t112;
  double t114;
  double t142;
  double t143;
  double t145;
  double t146;
  double t147;
  double t53;
  double t55;
  double t56;
  double t62;
  double t63;
  double t64;
  double t72;
  double t162;
  double t163;
  double t113;
  double t115;
  double t116;
  double t165;
  double t166;
  double t167;
  double t118;
  double t119;
  double t120;
  double t126;
  double t184;
  double t185;
  double t186;
  double t187;
  double t188;
  double t189;
  double t190;
  double t191;
  double t192;
  double t193;
  double t153;
  double t211;
  double t212;
  double t66;
  double t203;
  double t200;
  double t222;
  double t223;
  double t229;
  double t230;
  double t232;
  double t233;
  double t234;
  double t68;
  double t69;
  double t100;
  double t104;
  double t250;
  double t251;
  double t252;
  double t253;
  double t254;
  double t255;
  double t256;
  double t117;
  double t121;
  double t122;
  double t123;
  double t124;
  double t125;
  double t127;
  double t132;
  double t133;
  double t164;
  double t168;
  double t263;
  double t264;
  double t171;
  double t172;
  double t173;
  double t174;
  double t175;
  double t176;
  double t177;
  double t178;
  double t179;
  double t271;
  double t272;
  double t273;
  double t274;
  double t275;
  double t276;
  double t278;
  double t279;
  t9 = Sin(var1[2]);
  t6 = Cos(var1[2]);
  t20 = Cos(var1[3]);
  t21 = -1.*t20;
  t22 = 1. + t21;
  t24 = Sin(var1[3]);
  t23 = -0.0265*t22;
  t25 = -0.0695*t24;
  t28 = t23 + t25;
  t29 = -1.*t9*t28;
  t30 = -0.0695*t22;
  t32 = 0.0265*t24;
  t33 = t30 + t32;
  t34 = t6*t33;
  t38 = -1.*t20*t9;
  t40 = t6*t24;
  t41 = t38 + t40;
  t43 = t6*t20;
  t44 = t9*t24;
  t45 = t43 + t44;
  t50 = Cos(var1[4]);
  t51 = -1.*t50;
  t52 = 1. + t51;
  t54 = Sin(var1[4]);
  t84 = Cos(var1[5]);
  t85 = -1.*t84;
  t86 = 1. + t85;
  t88 = Sin(var1[5]);
  t87 = -0.0695*t86;
  t89 = -0.0265*t88;
  t90 = t87 + t89;
  t92 = t6*t90;
  t93 = -0.0265*t86;
  t94 = 0.0695*t88;
  t95 = t93 + t94;
  t96 = -1.*t9*t95;
  t101 = t6*t84;
  t102 = -1.*t9*t88;
  t103 = t101 + t102;
  t97 = -1.*t84*t9;
  t98 = -1.*t6*t88;
  t99 = t97 + t98;
  t108 = Cos(var1[6]);
  t110 = -1.*t108;
  t112 = 1. + t110;
  t114 = Sin(var1[6]);
  t142 = -1.*t6*t28;
  t143 = -1.*t9*t33;
  t145 = -1.*t6*t20;
  t146 = -1.*t9*t24;
  t147 = t145 + t146;
  t53 = -0.0265*t52;
  t55 = -0.2375*t54;
  t56 = t53 + t55;
  t62 = -0.2375*t52;
  t63 = 0.0265*t54;
  t64 = t62 + t63;
  t72 = t50*t41;
  t162 = -1.*t9*t90;
  t163 = -1.*t6*t95;
  t113 = -0.2375*t112;
  t115 = -0.0265*t114;
  t116 = t113 + t115;
  t165 = -1.*t6*t84;
  t166 = t9*t88;
  t167 = t165 + t166;
  t118 = -0.0265*t112;
  t119 = 0.2375*t114;
  t120 = t118 + t119;
  t126 = t108*t99;
  t184 = 0.0265*t20;
  t185 = t184 + t25;
  t186 = t9*t185;
  t187 = -0.0695*t20;
  t188 = -0.0265*t24;
  t189 = t187 + t188;
  t190 = t6*t189;
  t191 = t20*t9;
  t192 = -1.*t6*t24;
  t193 = t191 + t192;
  t153 = t50*t147;
  t211 = t6*t185;
  t212 = -1.*t9*t189;
  t66 = t50*t45;
  t203 = t50*t193;
  t200 = -1.*t193*t54;
  t222 = -1.*t45*t54;
  t223 = t203 + t222;
  t229 = 0.0265*t50;
  t230 = t229 + t55;
  t232 = -0.2375*t50;
  t233 = -0.0265*t54;
  t234 = t232 + t233;
  t68 = -1.*t41*t54;
  t69 = t66 + t68;
  t100 = -0.025413*t99;
  t104 = -0.15232*t103;
  t250 = -0.0265*t84;
  t251 = -0.0695*t88;
  t252 = t250 + t251;
  t253 = t9*t252;
  t254 = 0.0695*t84;
  t255 = t254 + t89;
  t256 = t6*t255;
  t117 = t103*t116;
  t121 = t99*t120;
  t122 = t108*t103;
  t123 = t99*t114;
  t124 = t122 + t123;
  t125 = -0.314506*t124;
  t127 = -1.*t103*t114;
  t132 = t126 + t127;
  t133 = -0.025226*t132;
  t164 = -0.15232*t99;
  t168 = -0.025413*t167;
  t263 = t6*t252;
  t264 = -1.*t9*t255;
  t171 = t99*t116;
  t172 = t167*t120;
  t173 = t108*t167;
  t174 = -1.*t99*t114;
  t175 = t173 + t174;
  t176 = -0.025226*t175;
  t177 = t167*t114;
  t178 = t126 + t177;
  t179 = -0.314506*t178;
  t271 = t84*t9;
  t272 = t6*t88;
  t273 = t271 + t272;
  t274 = -0.0265*t108;
  t275 = -0.2375*t114;
  t276 = t274 + t275;
  t278 = 0.2375*t108;
  t279 = t278 + t115;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=1.;
  p_output1[6]=0.2996793431028799*(0.69051*(t29 + t34 - 0.025367*t41 - 0.15232*t45) + 0.19964*(t29 + t34 + t41*t56 + t45*t64 - 0.314514*t69 - 0.025229*(t45*t54 + t72)) + 1.5566*(-0.046589*t6 + 0.026461*t9) + 0.69051*(t100 + t104 + t92 + t96) + 0.19964*(t117 + t121 + t125 + t133 + t92 + t96));
  p_output1[7]=0;
  p_output1[8]=0.2996793431028799*(0.69051*(t162 + t163 + t164 + t168) + 0.19964*(t162 + t163 + t171 + t172 + t176 + t179) + 0.69051*(t142 + t143 - 0.025367*t147 - 0.15232*t41) + 0.19964*(t142 + t143 - 0.025229*(t153 + t41*t54) + t147*t56 + t41*t64 - 0.314514*(-1.*t147*t54 + t72)) + 1.5566*(0.026461*t6 + 0.046589*t9));
  p_output1[9]=0.2996793431028799*(0.69051*(-0.15232*t147 + t186 + t190 - 0.025367*t193) + 0.19964*(t186 + t190 - 0.314514*(t153 + t200) - 0.025229*(t203 + t147*t54) + t193*t56 + t147*t64));
  p_output1[10]=0;
  p_output1[11]=0.2996793431028799*(0.69051*(-0.15232*t193 + t211 + t212 - 0.025367*t45) + 0.19964*(t211 + t212 - 0.314514*t223 + t45*t56 + t193*t64 - 0.025229*(t193*t54 + t66)));
  p_output1[12]=0.05982798405705895*(-0.025229*t223 + t193*t230 + t234*t45 - 0.314514*(t200 - 1.*t45*t50));
  p_output1[13]=0;
  p_output1[14]=0.05982798405705895*(t234*t41 + t230*t45 - 0.314514*(t222 - 1.*t41*t50) - 0.025229*t69);
  p_output1[15]=0.2996793431028799*(0.69051*(t100 + t104 + t253 + t256) + 0.19964*(t117 + t121 + t125 + t133 + t253 + t256));
  p_output1[16]=0;
  p_output1[17]=0.2996793431028799*(0.69051*(t164 + t168 + t263 + t264) + 0.19964*(t171 + t172 + t176 + t179 + t263 + t264));
  p_output1[18]=0.05982798405705895*(-0.025226*(t127 - 1.*t108*t273) - 0.314506*(t122 - 1.*t114*t273) + t273*t276 + t103*t279);
  p_output1[19]=0;
  p_output1[20]=0.05982798405705895*(-0.314506*t132 - 0.025226*(-1.*t103*t108 + t174) + t103*t276 + t279*t99);
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 3, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "J_COMPosition.hh"

namespace SymFunction
{

void J_COMPosition_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
