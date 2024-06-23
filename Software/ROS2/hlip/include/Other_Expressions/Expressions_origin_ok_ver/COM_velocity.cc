/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:24 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t18;
  double t8;
  double t10;
  double t23;
  double t32;
  double t40;
  double t46;
  double t48;
  double t49;
  double t16;
  double t25;
  double t30;
  double t78;
  double t82;
  double t87;
  double t89;
  double t92;
  double t93;
  double t94;
  double t100;
  double t42;
  double t104;
  double t105;
  double t106;
  double t117;
  double t118;
  double t63;
  double t66;
  double t144;
  double t149;
  double t152;
  double t154;
  double t158;
  double t159;
  double t160;
  double t148;
  double t150;
  double t151;
  double t179;
  double t180;
  double t181;
  double t182;
  double t183;
  double t184;
  double t185;
  double t186;
  double t164;
  double t187;
  double t188;
  double t190;
  double t197;
  double t198;
  double t168;
  double t173;
  double t223;
  double t224;
  double t225;
  double t226;
  double t227;
  double t228;
  double t229;
  double t230;
  double t231;
  double t232;
  double t233;
  double t234;
  double t119;
  double t120;
  double t123;
  double t124;
  double t125;
  double t251;
  double t252;
  double t191;
  double t194;
  double t253;
  double t254;
  double t256;
  double t257;
  double t258;
  double t259;
  double t260;
  double t199;
  double t201;
  double t202;
  double t203;
  double t204;
  double t205;
  double t206;
  double t207;
  double t208;
  double t209;
  double t210;
  double t211;
  double t213;
  double t34;
  double t43;
  double t53;
  double t55;
  double t56;
  double t241;
  double t242;
  double t243;
  double t69;
  double t278;
  double t279;
  double t70;
  double t153;
  double t155;
  double t156;
  double t161;
  double t166;
  double t303;
  double t304;
  double t306;
  double t307;
  double t308;
  double t297;
  double t328;
  double t329;
  double t132;
  double t245;
  double t305;
  double t309;
  double t344;
  double t345;
  double t312;
  double t313;
  double t314;
  double t315;
  double t316;
  double t317;
  double t318;
  double t319;
  t18 = Cos(var1[2]);
  t8 = Cos(var1[3]);
  t10 = Sin(var1[2]);
  t23 = Sin(var1[3]);
  t32 = Cos(var1[4]);
  t40 = Sin(var1[4]);
  t46 = t18*t8;
  t48 = t10*t23;
  t49 = t46 + t48;
  t16 = t8*t10;
  t25 = -1.*t18*t23;
  t30 = t16 + t25;
  t78 = 0.0265*t8;
  t82 = -0.0695*t23;
  t87 = t78 + t82;
  t89 = t10*t87;
  t92 = -0.0695*t8;
  t93 = -0.0265*t23;
  t94 = t92 + t93;
  t100 = t18*t94;
  t42 = -0.2375*t40;
  t104 = -1.*t18*t8;
  t105 = -1.*t10*t23;
  t106 = t104 + t105;
  t117 = -1.*t32;
  t118 = 1. + t117;
  t63 = -1.*t30*t40;
  t66 = t32*t30;
  t144 = Cos(var1[5]);
  t149 = Sin(var1[5]);
  t152 = Cos(var1[6]);
  t154 = Sin(var1[6]);
  t158 = t18*t144;
  t159 = -1.*t10*t149;
  t160 = t158 + t159;
  t148 = t144*t10;
  t150 = t18*t149;
  t151 = t148 + t150;
  t179 = -0.0265*t144;
  t180 = -0.0695*t149;
  t181 = t179 + t180;
  t182 = t10*t181;
  t183 = 0.0695*t144;
  t184 = -0.0265*t149;
  t185 = t183 + t184;
  t186 = t18*t185;
  t164 = -0.0265*t154;
  t187 = -1.*t144*t10;
  t188 = -1.*t18*t149;
  t190 = t187 + t188;
  t197 = -1.*t152;
  t198 = 1. + t197;
  t168 = t152*t160;
  t173 = -1.*t160*t154;
  t223 = -1.*t8;
  t224 = 1. + t223;
  t225 = -0.0265*t224;
  t226 = t225 + t82;
  t227 = -1.*t10*t226;
  t228 = -0.0695*t224;
  t229 = 0.0265*t23;
  t230 = t228 + t229;
  t231 = t18*t230;
  t232 = -1.*t8*t10;
  t233 = t18*t23;
  t234 = t232 + t233;
  t119 = -0.0265*t118;
  t120 = t119 + t42;
  t123 = -0.2375*t118;
  t124 = 0.0265*t40;
  t125 = t123 + t124;
  t251 = -1.*t144;
  t252 = 1. + t251;
  t191 = -0.025413*t190;
  t194 = -0.15232*t160;
  t253 = -0.0695*t252;
  t254 = t253 + t184;
  t256 = t18*t254;
  t257 = -0.0265*t252;
  t258 = 0.0695*t149;
  t259 = t257 + t258;
  t260 = -1.*t10*t259;
  t199 = -0.2375*t198;
  t201 = t199 + t164;
  t202 = t160*t201;
  t203 = -0.0265*t198;
  t204 = 0.2375*t154;
  t205 = t203 + t204;
  t206 = t190*t205;
  t207 = t190*t154;
  t208 = t168 + t207;
  t209 = -0.314506*t208;
  t210 = t152*t190;
  t211 = t210 + t173;
  t213 = -0.025226*t211;
  t34 = 0.0265*t32;
  t43 = t34 + t42;
  t53 = -0.2375*t32;
  t55 = -0.0265*t40;
  t56 = t53 + t55;
  t241 = t32*t49;
  t242 = -1.*t234*t40;
  t243 = t241 + t242;
  t69 = -1.*t49*t40;
  t278 = t18*t87;
  t279 = -1.*t10*t94;
  t70 = t66 + t69;
  t153 = -0.0265*t152;
  t155 = -0.2375*t154;
  t156 = t153 + t155;
  t161 = 0.2375*t152;
  t166 = t161 + t164;
  t303 = t18*t181;
  t304 = -1.*t10*t185;
  t306 = -1.*t18*t144;
  t307 = t10*t149;
  t308 = t306 + t307;
  t297 = -1.*t190*t154;
  t328 = -1.*t18*t226;
  t329 = -1.*t10*t230;
  t132 = t32*t106;
  t245 = t32*t234;
  t305 = -0.15232*t190;
  t309 = -0.025413*t308;
  t344 = -1.*t10*t254;
  t345 = -1.*t18*t259;
  t312 = t190*t201;
  t313 = t308*t205;
  t314 = t152*t308;
  t315 = t314 + t297;
  t316 = -0.025226*t315;
  t317 = t308*t154;
  t318 = t210 + t317;
  t319 = -0.314506*t318;
  p_output1[0]=var2[0] + 0.2996793431028799*(1.5566*(0.026461*t10 - 0.046589*t18) + 0.69051*(t191 + t194 + t256 + t260) + 0.19964*(t202 + t206 + t209 + t213 + t256 + t260) + 0.69051*(t227 + t231 - 0.025367*t234 - 0.15232*t49) + 0.19964*(t227 + t231 + t120*t234 - 0.314514*t243 + t125*t49 - 0.025229*(t245 + t40*t49)))*var2[2] + 0.2996793431028799*(0.69051*(t100 - 0.15232*t106 - 0.025367*t30 + t89) + 0.19964*(t100 + t106*t125 + t120*t30 - 0.314514*(t132 + t63) - 0.025229*(t106*t40 + t66) + t89))*var2[3] + 0.05982798405705895*(t30*t43 + t49*t56 - 0.314514*(-1.*t32*t49 + t63) - 0.025229*t70)*var2[4] + 0.2996793431028799*(0.69051*(t182 + t186 + t191 + t194) + 0.19964*(t182 + t186 + t202 + t206 + t209 + t213))*var2[5] + 0.05982798405705895*(t151*t156 + t160*t166 - 0.314506*(-1.*t151*t154 + t168) - 0.025226*(-1.*t151*t152 + t173))*var2[6];
  p_output1[1]=0;
  p_output1[2]=var2[1] + 0.2996793431028799*(1.5566*(0.046589*t10 + 0.026461*t18) + 0.69051*(-0.025367*t106 - 0.15232*t234 + t328 + t329) + 0.69051*(t305 + t309 + t344 + t345) + 0.19964*(t312 + t313 + t316 + t319 + t344 + t345) + 0.19964*(t106*t120 + t125*t234 + t328 + t329 - 0.314514*(t245 - 1.*t106*t40) - 0.025229*(t132 + t234*t40)))*var2[2] + 0.2996793431028799*(0.69051*(t278 + t279 - 0.15232*t30 - 0.025367*t49) + 0.19964*(t278 + t279 + t125*t30 - 0.025229*(t241 + t30*t40) + t120*t49 - 0.314514*t70))*var2[3] + 0.05982798405705895*(-0.025229*t243 + t43*t49 + t234*t56 - 0.314514*(-1.*t234*t32 + t69))*var2[4] + 0.2996793431028799*(0.69051*(t303 + t304 + t305 + t309) + 0.19964*(t303 + t304 + t312 + t313 + t316 + t319))*var2[5] + 0.05982798405705895*(t156*t160 + t166*t190 - 0.314506*t211 - 0.025226*(-1.*t152*t160 + t297))*var2[6];
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

  double *var1,*var2;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 2)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "Two input(s) required (var1,var2).");
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
  mrows = mxGetM(prhs[1]);
  ncols = mxGetN(prhs[1]);
  if( !mxIsDouble(prhs[1]) || mxIsComplex(prhs[1]) ||
    ( !(mrows == 7 && ncols == 1) && 
      !(mrows == 1 && ncols == 7))) 
    {
      mexErrMsgIdAndTxt( "MATLAB:MShaped:inputNotRealVector", "var2 is wrong.");
    }

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
  var2 = mxGetPr(prhs[1]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 3, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "COM_velocity.hh"

namespace SymFunction
{

void COM_velocity_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
