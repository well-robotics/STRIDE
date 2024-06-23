/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:36:33 GMT-05:00
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
  double t6;
  double t3;
  double t9;
  double t11;
  double t21;
  double t23;
  double t27;
  double t28;
  double t29;
  double t42;
  double t44;
  double t54;
  double t56;
  double t60;
  double t61;
  double t62;
  double t16;
  double t17;
  double t18;
  double t10;
  double t12;
  double t13;
  double t22;
  double t24;
  double t25;
  double t26;
  double t30;
  double t31;
  double t34;
  double t35;
  double t36;
  double t37;
  double t38;
  double t39;
  double t49;
  double t50;
  double t51;
  double t43;
  double t45;
  double t46;
  double t55;
  double t57;
  double t58;
  double t59;
  double t63;
  double t64;
  double t67;
  double t68;
  double t69;
  double t70;
  double t71;
  double t72;
  double t107;
  double t108;
  double t109;
  double t103;
  double t104;
  double t105;
  double t117;
  double t118;
  double t119;
  double t121;
  double t122;
  double t123;
  double t124;
  double t125;
  double t126;
  double t130;
  double t131;
  double t132;
  double t120;
  double t127;
  double t128;
  double t144;
  double t145;
  double t146;
  double t140;
  double t141;
  double t142;
  double t154;
  double t155;
  double t156;
  double t158;
  double t159;
  double t160;
  double t161;
  double t162;
  double t163;
  double t167;
  double t168;
  double t169;
  double t157;
  double t164;
  double t165;
  double t76;
  double t77;
  double t78;
  double t79;
  double t80;
  double t81;
  double t82;
  double t83;
  double t84;
  double t85;
  double t86;
  double t87;
  double t88;
  double t89;
  double t90;
  double t91;
  double t92;
  double t93;
  double t94;
  double t95;
  double t96;
  double t97;
  double t98;
  double t99;
  double t100;
  double t4;
  double t5;
  double t7;
  double t8;
  double t19;
  double t20;
  double t52;
  double t53;
  double t106;
  double t110;
  double t111;
  double t113;
  double t114;
  double t115;
  double t129;
  double t133;
  double t134;
  double t136;
  double t137;
  double t138;
  double t143;
  double t147;
  double t148;
  double t150;
  double t151;
  double t152;
  double t166;
  double t170;
  double t171;
  double t173;
  double t174;
  double t175;
  double t184;
  double t185;
  double t186;
  double t180;
  double t181;
  double t182;
  double t198;
  double t199;
  double t200;
  double t194;
  double t195;
  double t196;
  double t101;
  double t102;
  double t112;
  double t116;
  double t135;
  double t139;
  double t149;
  double t153;
  double t172;
  double t176;
  double t177;
  double t219;
  double t220;
  double t221;
  double t222;
  double t223;
  double t224;
  double t225;
  double t226;
  double t227;
  double t228;
  double t229;
  double t178;
  double t179;
  double t183;
  double t187;
  double t188;
  double t230;
  double t231;
  double t232;
  double t233;
  double t234;
  double t263;
  double t264;
  double t265;
  double t266;
  double t267;
  double t189;
  double t190;
  double t191;
  double t235;
  double t236;
  double t237;
  double t268;
  double t269;
  double t270;
  double t284;
  double t285;
  double t286;
  double t192;
  double t193;
  double t197;
  double t201;
  double t202;
  double t238;
  double t239;
  double t240;
  double t241;
  double t242;
  double t271;
  double t272;
  double t273;
  double t274;
  double t275;
  double t203;
  double t204;
  double t205;
  double t243;
  double t244;
  double t245;
  double t276;
  double t277;
  double t278;
  double t292;
  double t293;
  double t294;
  t6 = Sin(var1[2]);
  t3 = Cos(var1[2]);
  t9 = Cos(var1[3]);
  t11 = Sin(var1[3]);
  t21 = Cos(var1[4]);
  t23 = Sin(var1[4]);
  t27 = t9*t21;
  t28 = -1.*t11*t23;
  t29 = t27 + t28;
  t42 = Cos(var1[5]);
  t44 = Sin(var1[5]);
  t54 = Cos(var1[6]);
  t56 = Sin(var1[6]);
  t60 = t42*t54;
  t61 = -1.*t44*t56;
  t62 = t60 + t61;
  t16 = t3*t9;
  t17 = -1.*t6*t11;
  t18 = t16 + t17;
  t10 = t9*t6;
  t12 = t3*t11;
  t13 = t10 + t12;
  t22 = -1.*t21*t11;
  t24 = -1.*t9*t23;
  t25 = t22 + t24;
  t26 = t6*t25;
  t30 = t3*t29;
  t31 = t26 + t30;
  t34 = t21*t11;
  t35 = t9*t23;
  t36 = t34 + t35;
  t37 = t3*t36;
  t38 = t6*t29;
  t39 = t37 + t38;
  t49 = t3*t42;
  t50 = -1.*t6*t44;
  t51 = t49 + t50;
  t43 = t42*t6;
  t45 = t3*t44;
  t46 = t43 + t45;
  t55 = -1.*t54*t44;
  t57 = -1.*t42*t56;
  t58 = t55 + t57;
  t59 = t6*t58;
  t63 = t3*t62;
  t64 = t59 + t63;
  t67 = t54*t44;
  t68 = t42*t56;
  t69 = t67 + t68;
  t70 = t3*t69;
  t71 = t6*t62;
  t72 = t70 + t71;
  t107 = -0.00102*t9;
  t108 = -0.078722*t11;
  t109 = t107 + t108;
  t103 = -0.078722*t9;
  t104 = 0.00102*t11;
  t105 = t103 + t104;
  t117 = -0.022663*t21;
  t118 = -0.007370999999999989*t23;
  t119 = t117 + t118;
  t121 = -1.*t21;
  t122 = 1. + t121;
  t123 = -0.16*t122;
  t124 = -0.167371*t21;
  t125 = 0.022663*t23;
  t126 = t123 + t124 + t125;
  t130 = -1.*t11*t119;
  t131 = t9*t126;
  t132 = t130 + t131;
  t120 = t9*t119;
  t127 = t11*t126;
  t128 = t120 + t127;
  t144 = -0.001112*t42;
  t145 = -0.078865*t44;
  t146 = t144 + t145;
  t140 = -0.078865*t42;
  t141 = 0.001112*t44;
  t142 = t140 + t141;
  t154 = -0.022659*t54;
  t155 = -0.007367999999999986*t56;
  t156 = t154 + t155;
  t158 = -1.*t54;
  t159 = 1. + t158;
  t160 = -0.16*t159;
  t161 = -0.167368*t54;
  t162 = 0.022659*t56;
  t163 = t160 + t161 + t162;
  t167 = -1.*t44*t156;
  t168 = t42*t163;
  t169 = t167 + t168;
  t157 = t42*t156;
  t164 = t44*t163;
  t165 = t157 + t164;
  t76 = -1.*t9*t6;
  t77 = -1.*t3*t11;
  t78 = t76 + t77;
  t79 = 0.64654*t78*t18;
  t80 = 0.64654*t13*t18;
  t81 = t3*t25;
  t82 = -1.*t6*t29;
  t83 = t81 + t82;
  t84 = 0.14994*t31*t83;
  t85 = -1.*t6*t36;
  t86 = t85 + t30;
  t87 = 0.14994*t86*t39;
  t88 = -1.*t42*t6;
  t89 = -1.*t3*t44;
  t90 = t88 + t89;
  t91 = 0.64429*t90*t51;
  t92 = 0.64429*t46*t51;
  t93 = t3*t58;
  t94 = -1.*t6*t62;
  t95 = t93 + t94;
  t96 = 0.14994*t64*t95;
  t97 = -1.*t6*t69;
  t98 = t97 + t63;
  t99 = 0.14994*t98*t72;
  t100 = 0. + t79 + t80 + t84 + t87 + t91 + t92 + t96 + t99;
  t4 = Power(t3,2);
  t5 = 1.2497*t4;
  t7 = Power(t6,2);
  t8 = 1.2497*t7;
  t19 = Power(t18,2);
  t20 = 0.64654*t19;
  t52 = Power(t51,2);
  t53 = 0.64429*t52;
  t106 = t9*t105;
  t110 = t109*t11;
  t111 = t106 + t110;
  t113 = -1.*t9*t109;
  t114 = t105*t11;
  t115 = t113 + t114;
  t129 = t128*t36;
  t133 = t132*t29;
  t134 = t129 + t133;
  t136 = -1.*t132*t25;
  t137 = -1.*t128*t29;
  t138 = t136 + t137;
  t143 = t42*t142;
  t147 = t146*t44;
  t148 = t143 + t147;
  t150 = -1.*t42*t146;
  t151 = t142*t44;
  t152 = t150 + t151;
  t166 = t165*t69;
  t170 = t169*t62;
  t171 = t166 + t170;
  t173 = -1.*t169*t58;
  t174 = -1.*t165*t62;
  t175 = t173 + t174;
  t184 = -1.*t21*t119;
  t185 = t126*t23;
  t186 = t184 + t185;
  t180 = t21*t126;
  t181 = t119*t23;
  t182 = t180 + t181;
  t198 = -1.*t54*t156;
  t199 = t163*t56;
  t200 = t198 + t199;
  t194 = t54*t163;
  t195 = t156*t56;
  t196 = t194 + t195;
  t101 = 0.0099626084*t3;
  t102 = -0.0000524874*t6;
  t112 = 0.64654*t18*t111;
  t116 = 0.64654*t13*t115;
  t135 = 0.14994*t31*t134;
  t139 = 0.14994*t39*t138;
  t149 = 0.64429*t51*t148;
  t153 = 0.64429*t46*t152;
  t172 = 0.14994*t64*t171;
  t176 = 0.14994*t72*t175;
  t177 = t101 + t102 + t112 + t116 + t135 + t139 + t149 + t153 + t172 + t176;
  t219 = -0.0000524874*t3;
  t220 = -0.0099626084*t6;
  t221 = 0.64654*t78*t111;
  t222 = 0.64654*t18*t115;
  t223 = 0.14994*t83*t134;
  t224 = 0.14994*t86*t138;
  t225 = 0.64429*t90*t148;
  t226 = 0.64429*t51*t152;
  t227 = 0.14994*t95*t171;
  t228 = 0.14994*t98*t175;
  t229 = t219 + t220 + t221 + t222 + t223 + t224 + t225 + t226 + t227 + t228;
  t178 = 0.0006594708000000001*t13;
  t179 = -0.05089692188*t18;
  t183 = 0.14994*t182*t31;
  t187 = 0.14994*t186*t39;
  t188 = t178 + t179 + t183 + t187;
  t230 = -0.05089692188*t78;
  t231 = 0.0006594708000000001*t18;
  t232 = 0.14994*t186*t86;
  t233 = 0.14994*t182*t83;
  t234 = t230 + t231 + t232 + t233;
  t263 = -0.05089692188*t111;
  t264 = 0.0006594708000000001*t115;
  t265 = 0.14994*t182*t134;
  t266 = 0.14994*t186*t138;
  t267 = 0.000225 + t263 + t264 + t265 + t266;
  t189 = -0.0011052077399999983*t31;
  t190 = 0.0033980902199999994*t39;
  t191 = t189 + t190;
  t235 = 0.0033980902199999994*t86;
  t236 = -0.0011052077399999983*t83;
  t237 = t235 + t236;
  t268 = -0.0011052077399999983*t134;
  t269 = 0.0033980902199999994*t138;
  t270 = 0.000088 + t268 + t269;
  t284 = -0.0011052077399999983*t182;
  t285 = 0.0033980902199999994*t186;
  t286 = 0.000088 + t284 + t285;
  t192 = 0.00071645048*t46;
  t193 = -0.050811930850000006*t51;
  t197 = 0.14994*t196*t64;
  t201 = 0.14994*t200*t72;
  t202 = t192 + t193 + t197 + t201;
  t238 = -0.050811930850000006*t90;
  t239 = 0.00071645048*t51;
  t240 = 0.14994*t200*t98;
  t241 = 0.14994*t196*t95;
  t242 = t238 + t239 + t240 + t241;
  t271 = -0.050811930850000006*t148;
  t272 = 0.00071645048*t152;
  t273 = 0.14994*t196*t171;
  t274 = 0.14994*t200*t175;
  t275 = 0.000225 + t271 + t272 + t273 + t274;
  t203 = -0.0011047579199999977*t64;
  t204 = 0.0033974904599999994*t72;
  t205 = t203 + t204;
  t243 = 0.0033974904599999994*t98;
  t244 = -0.0011047579199999977*t95;
  t245 = t243 + t244;
  t276 = -0.0011047579199999977*t171;
  t277 = 0.0033974904599999994*t175;
  t278 = 0.000088 + t276 + t277;
  t292 = -0.0011047579199999977*t196;
  t293 = 0.0033974904599999994*t200;
  t294 = 0.000088 + t292 + t293;
  p_output1[0]=0.64654*Power(t13,2) + t20 + 0.14994*Power(t31,2) + 0.14994*Power(t39,2) + 0.64429*Power(t46,2) + t5 + t53 + 0.14994*Power(t64,2) + 0.14994*Power(t72,2) + t8;
  p_output1[1]=t100;
  p_output1[2]=t177;
  p_output1[3]=t188;
  p_output1[4]=t191;
  p_output1[5]=t202;
  p_output1[6]=t205;
  p_output1[7]=t100;
  p_output1[8]=t20 + t5 + t53 + 0.64654*Power(t78,2) + t8 + 0.14994*Power(t83,2) + 0.14994*Power(t86,2) + 0.64429*Power(t90,2) + 0.14994*Power(t95,2) + 0.14994*Power(t98,2);
  p_output1[9]=t229;
  p_output1[10]=t234;
  p_output1[11]=t237;
  p_output1[12]=t242;
  p_output1[13]=t245;
  p_output1[14]=t177;
  p_output1[15]=t229;
  p_output1[16]=0.0009224241186355999 + 0.64654*Power(t111,2) + 0.64654*Power(t115,2) + 0.14994*Power(t134,2) + 0.14994*Power(t138,2) + 0.64429*Power(t148,2) + 0.64429*Power(t152,2) + 0.14994*Power(t171,2) + 0.14994*Power(t175,2);
  p_output1[17]=t267;
  p_output1[18]=t270;
  p_output1[19]=t275;
  p_output1[20]=t278;
  p_output1[21]=t188;
  p_output1[22]=t234;
  p_output1[23]=t267;
  p_output1[24]=0.8342323801444533 + 0.14994*Power(t182,2) + 0.14994*Power(t186,2);
  p_output1[25]=t286;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=t191;
  p_output1[29]=t237;
  p_output1[30]=t270;
  p_output1[31]=t286;
  p_output1[32]=0.8301731574049074;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=t202;
  p_output1[36]=t242;
  p_output1[37]=t275;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.834233079619419 + 0.14994*Power(t196,2) + 0.14994*Power(t200,2);
  p_output1[41]=t294;
  p_output1[42]=t205;
  p_output1[43]=t245;
  p_output1[44]=t278;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=t294;
  p_output1[48]=0.8301731235926877;
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

#include "Mmat_five_link_walker.hh"

namespace SymFunction
{

void Mmat_five_link_walker_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
