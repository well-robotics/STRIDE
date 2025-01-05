/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:28:59 GMT-05:00
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
  double t10;
  double t12;
  double t13;
  double t16;
  double t17;
  double t18;
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
  double t103;
  double t104;
  double t110;
  double t111;
  double t112;
  double t113;
  double t105;
  double t106;
  double t107;
  double t108;
  double t122;
  double t123;
  double t124;
  double t125;
  double t126;
  double t127;
  double t129;
  double t130;
  double t131;
  double t132;
  double t136;
  double t137;
  double t138;
  double t139;
  double t121;
  double t128;
  double t133;
  double t134;
  double t147;
  double t148;
  double t154;
  double t155;
  double t156;
  double t157;
  double t149;
  double t150;
  double t151;
  double t152;
  double t166;
  double t167;
  double t168;
  double t169;
  double t170;
  double t171;
  double t173;
  double t174;
  double t175;
  double t176;
  double t180;
  double t181;
  double t182;
  double t183;
  double t165;
  double t172;
  double t177;
  double t178;
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
  double t109;
  double t114;
  double t115;
  double t117;
  double t118;
  double t119;
  double t135;
  double t140;
  double t141;
  double t143;
  double t144;
  double t145;
  double t153;
  double t158;
  double t159;
  double t161;
  double t162;
  double t163;
  double t179;
  double t184;
  double t185;
  double t187;
  double t188;
  double t189;
  double t200;
  double t201;
  double t202;
  double t203;
  double t204;
  double t194;
  double t195;
  double t196;
  double t197;
  double t198;
  double t218;
  double t219;
  double t220;
  double t221;
  double t222;
  double t212;
  double t213;
  double t214;
  double t215;
  double t216;
  double t101;
  double t102;
  double t116;
  double t120;
  double t142;
  double t146;
  double t160;
  double t164;
  double t186;
  double t190;
  double t191;
  double t241;
  double t242;
  double t243;
  double t244;
  double t245;
  double t246;
  double t247;
  double t248;
  double t249;
  double t250;
  double t251;
  double t192;
  double t193;
  double t199;
  double t205;
  double t206;
  double t252;
  double t253;
  double t254;
  double t255;
  double t256;
  double t285;
  double t286;
  double t287;
  double t288;
  double t289;
  double t207;
  double t208;
  double t209;
  double t257;
  double t258;
  double t259;
  double t290;
  double t291;
  double t292;
  double t306;
  double t307;
  double t308;
  double t210;
  double t211;
  double t217;
  double t223;
  double t224;
  double t260;
  double t261;
  double t262;
  double t263;
  double t264;
  double t293;
  double t294;
  double t295;
  double t296;
  double t297;
  double t225;
  double t226;
  double t227;
  double t265;
  double t266;
  double t267;
  double t298;
  double t299;
  double t300;
  double t314;
  double t315;
  double t316;
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
  t10 = t9*t6;
  t12 = -1.*t3*t11;
  t13 = t10 + t12;
  t16 = t3*t9;
  t17 = t6*t11;
  t18 = t16 + t17;
  t22 = t21*t11;
  t24 = t9*t23;
  t25 = t22 + t24;
  t26 = t6*t25;
  t30 = t3*t29;
  t31 = t26 + t30;
  t34 = -1.*t21*t11;
  t35 = -1.*t9*t23;
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
  t103 = -1.*t9;
  t104 = 1. + t103;
  t110 = -0.0695*t104;
  t111 = -0.15232*t9;
  t112 = 0.0011329999999999986*t11;
  t113 = t110 + t111 + t112;
  t105 = -0.0265*t104;
  t106 = -0.025367*t9;
  t107 = 0.08282*t11;
  t108 = t105 + t106 + t107;
  t122 = -1.*t21;
  t123 = 1. + t122;
  t124 = -0.2375*t123;
  t125 = -0.314514*t21;
  t126 = 0.0012709999999999978*t23;
  t127 = t124 + t125 + t126;
  t129 = -0.0265*t123;
  t130 = -0.025229*t21;
  t131 = 0.07701400000000003*t23;
  t132 = t129 + t130 + t131;
  t136 = -0.0695*t11;
  t137 = -1.*t11*t127;
  t138 = t9*t132;
  t139 = t105 + t136 + t137 + t138;
  t121 = 0.0265*t11;
  t128 = t9*t127;
  t133 = t11*t132;
  t134 = t110 + t121 + t128 + t133;
  t147 = -1.*t42;
  t148 = 1. + t147;
  t154 = -0.0265*t148;
  t155 = -0.025413*t42;
  t156 = -0.08282*t44;
  t157 = t154 + t155 + t156;
  t149 = -0.0695*t148;
  t150 = -0.15232*t42;
  t151 = -0.0010869999999999977*t44;
  t152 = t149 + t150 + t151;
  t166 = -1.*t54;
  t167 = 1. + t166;
  t168 = -0.0265*t167;
  t169 = -0.025226*t54;
  t170 = -0.07700600000000002*t56;
  t171 = t168 + t169 + t170;
  t173 = -0.2375*t167;
  t174 = -0.314506*t54;
  t175 = -0.0012740000000000008*t56;
  t176 = t173 + t174 + t175;
  t180 = -0.0265*t44;
  t181 = -1.*t44*t171;
  t182 = t42*t176;
  t183 = t149 + t180 + t181 + t182;
  t165 = 0.0695*t44;
  t172 = t42*t171;
  t177 = t44*t176;
  t178 = t154 + t165 + t172 + t177;
  t76 = 0.69051*t13*t18;
  t77 = -1.*t9*t6;
  t78 = t3*t11;
  t79 = t77 + t78;
  t80 = 0.69051*t79*t18;
  t81 = t3*t25;
  t82 = -1.*t6*t29;
  t83 = t81 + t82;
  t84 = 0.19964*t31*t83;
  t85 = -1.*t6*t36;
  t86 = t85 + t30;
  t87 = 0.19964*t86*t39;
  t88 = -1.*t42*t6;
  t89 = -1.*t3*t44;
  t90 = t88 + t89;
  t91 = 0.69051*t90*t51;
  t92 = 0.69051*t46*t51;
  t93 = t3*t58;
  t94 = -1.*t6*t62;
  t95 = t93 + t94;
  t96 = 0.19964*t64*t95;
  t97 = -1.*t6*t69;
  t98 = t97 + t63;
  t99 = 0.19964*t98*t72;
  t100 = 0. + t76 + t80 + t84 + t87 + t91 + t92 + t96 + t99;
  t4 = Power(t3,2);
  t5 = 1.5566*t4;
  t7 = Power(t6,2);
  t8 = 1.5566*t7;
  t19 = Power(t18,2);
  t20 = 0.69051*t19;
  t52 = Power(t51,2);
  t53 = 0.69051*t52;
  t109 = -1.*t9*t108;
  t114 = -1.*t113*t11;
  t115 = t109 + t114;
  t117 = t9*t113;
  t118 = -1.*t108*t11;
  t119 = t117 + t118;
  t135 = -1.*t134*t25;
  t140 = -1.*t139*t29;
  t141 = t135 + t140;
  t143 = t139*t36;
  t144 = t134*t29;
  t145 = t143 + t144;
  t153 = t42*t152;
  t158 = t157*t44;
  t159 = t153 + t158;
  t161 = -1.*t42*t157;
  t162 = t152*t44;
  t163 = t161 + t162;
  t179 = t178*t69;
  t184 = t183*t62;
  t185 = t179 + t184;
  t187 = -1.*t183*t58;
  t188 = -1.*t178*t62;
  t189 = t187 + t188;
  t200 = 0.0265*t21;
  t201 = t21*t132;
  t202 = 0.0695*t23;
  t203 = t127*t23;
  t204 = t200 + t201 + t202 + t203;
  t194 = -0.0695*t21;
  t195 = -1.*t21*t127;
  t196 = 0.0265*t23;
  t197 = t132*t23;
  t198 = t194 + t195 + t196 + t197;
  t218 = -0.0265*t54;
  t219 = -1.*t54*t171;
  t220 = 0.0695*t56;
  t221 = t176*t56;
  t222 = t218 + t219 + t220 + t221;
  t212 = 0.0695*t54;
  t213 = t54*t176;
  t214 = 0.0265*t56;
  t215 = t171*t56;
  t216 = t212 + t213 + t214 + t215;
  t101 = -0.0725204374*t3;
  t102 = 0.0411891926*t6;
  t116 = 0.69051*t13*t115;
  t120 = 0.69051*t18*t119;
  t142 = 0.19964*t39*t141;
  t146 = 0.19964*t31*t145;
  t160 = 0.69051*t51*t159;
  t164 = 0.69051*t46*t163;
  t186 = 0.19964*t64*t185;
  t190 = 0.19964*t72*t189;
  t191 = t101 + t102 + t116 + t120 + t142 + t146 + t160 + t164 + t186 + t190;
  t241 = 0.0411891926*t3;
  t242 = 0.0725204374*t6;
  t243 = 0.69051*t18*t115;
  t244 = 0.69051*t79*t119;
  t245 = 0.19964*t86*t141;
  t246 = 0.19964*t83*t145;
  t247 = 0.69051*t90*t159;
  t248 = 0.69051*t51*t163;
  t249 = 0.19964*t95*t185;
  t250 = 0.19964*t98*t189;
  t251 = t241 + t242 + t243 + t244 + t245 + t246 + t247 + t248 + t249 + t250;
  t192 = 0.0007823478299999989*t13;
  t193 = 0.0571880382*t18;
  t199 = 0.19964*t198*t31;
  t205 = 0.19964*t204*t39;
  t206 = t192 + t193 + t199 + t205;
  t252 = 0.0571880382*t79;
  t253 = 0.0007823478299999989*t18;
  t254 = 0.19964*t204*t86;
  t255 = 0.19964*t198*t83;
  t256 = t252 + t253 + t254 + t255;
  t285 = 0.0007823478299999989*t115;
  t286 = 0.0571880382*t119;
  t287 = 0.19964*t204*t141;
  t288 = 0.19964*t198*t145;
  t289 = -0.000356 + t285 + t286 + t287 + t288;
  t207 = 0.015375074960000006*t31;
  t208 = 0.0002537424399999996*t39;
  t209 = t207 + t208;
  t257 = 0.0002537424399999996*t86;
  t258 = 0.015375074960000006*t83;
  t259 = t257 + t258;
  t290 = 0.0002537424399999996*t141;
  t291 = 0.015375074960000006*t145;
  t292 = -0.000106 + t290 + t291;
  t306 = 0.0002537424399999996*t204;
  t307 = 0.015375074960000006*t198;
  t308 = 0.000106 + t306 + t307;
  t210 = -0.0007505843699999984*t46;
  t211 = -0.0571880382*t51;
  t217 = 0.19964*t216*t64;
  t223 = 0.19964*t222*t72;
  t224 = t210 + t211 + t217 + t223;
  t260 = -0.0571880382*t90;
  t261 = -0.0007505843699999984*t51;
  t262 = 0.19964*t222*t98;
  t263 = 0.19964*t216*t95;
  t264 = t260 + t261 + t262 + t263;
  t293 = -0.0571880382*t159;
  t294 = -0.0007505843699999984*t163;
  t295 = 0.19964*t216*t185;
  t296 = 0.19964*t222*t189;
  t297 = 0.000355 + t293 + t294 + t295 + t296;
  t225 = -0.015373477840000005*t64;
  t226 = -0.0002543413600000002*t72;
  t227 = t225 + t226;
  t265 = -0.0002543413600000002*t98;
  t266 = -0.015373477840000005*t95;
  t267 = t265 + t266;
  t298 = -0.015373477840000005*t185;
  t299 = -0.0002543413600000002*t189;
  t300 = 0.000105 + t298 + t299;
  t314 = -0.015373477840000005*t216;
  t315 = -0.0002543413600000002*t222;
  t316 = 0.000105 + t314 + t315;
  p_output1[0]=0.69051*Power(t13,2) + t20 + 0.19964*Power(t31,2) + 0.19964*Power(t39,2) + 0.69051*Power(t46,2) + t5 + t53 + 0.19964*Power(t64,2) + 0.19964*Power(t72,2) + t8;
  p_output1[1]=t100;
  p_output1[2]=t191;
  p_output1[3]=t206;
  p_output1[4]=t209;
  p_output1[5]=t224;
  p_output1[6]=t227;
  p_output1[7]=t100;
  p_output1[8]=t20 + t5 + t53 + 0.69051*Power(t79,2) + t8 + 0.19964*Power(t83,2) + 0.19964*Power(t86,2) + 0.69051*Power(t90,2) + 0.19964*Power(t95,2) + 0.19964*Power(t98,2);
  p_output1[9]=t251;
  p_output1[10]=t256;
  p_output1[11]=t259;
  p_output1[12]=t264;
  p_output1[13]=t267;
  p_output1[14]=t191;
  p_output1[15]=t251;
  p_output1[16]=0.0058385618834172 + 0.69051*Power(t115,2) + 0.69051*Power(t119,2) + 0.19964*Power(t141,2) + 0.19964*Power(t145,2) + 0.69051*Power(t159,2) + 0.69051*Power(t163,2) + 0.19964*Power(t185,2) + 0.19964*Power(t189,2);
  p_output1[17]=t289;
  p_output1[18]=t292;
  p_output1[19]=t297;
  p_output1[20]=t300;
  p_output1[21]=t206;
  p_output1[22]=t256;
  p_output1[23]=t289;
  p_output1[24]=1.6881471997238156 + 0.19964*Power(t198,2) + 0.19964*Power(t204,2);
  p_output1[25]=t308;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=t209;
  p_output1[29]=t259;
  p_output1[30]=t292;
  p_output1[31]=t308;
  p_output1[32]=1.6843444185296108;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=t224;
  p_output1[36]=t264;
  p_output1[37]=t297;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=1.6881461292089341 + 0.19964*Power(t216,2) + 0.19964*Power(t222,2);
  p_output1[41]=t316;
  p_output1[42]=t227;
  p_output1[43]=t267;
  p_output1[44]=t300;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=t316;
  p_output1[48]=1.6843431740654398;
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

#include "Mmat_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Mmat_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
