/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:25 GMT-05:00
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
  double t34;
  double t43;
  double t16;
  double t44;
  double t53;
  double t56;
  double t42;
  double t46;
  double t48;
  double t65;
  double t66;
  double t70;
  double t119;
  double t120;
  double t121;
  double t123;
  double t124;
  double t125;
  double t126;
  double t132;
  double t60;
  double t133;
  double t134;
  double t135;
  double t140;
  double t141;
  double t82;
  double t100;
  double t177;
  double t179;
  double t182;
  double t184;
  double t178;
  double t180;
  double t181;
  double t195;
  double t196;
  double t199;
  double t216;
  double t217;
  double t218;
  double t219;
  double t220;
  double t221;
  double t225;
  double t228;
  double t202;
  double t231;
  double t233;
  double t235;
  double t239;
  double t240;
  double t206;
  double t209;
  double t271;
  double t272;
  double t274;
  double t275;
  double t276;
  double t277;
  double t278;
  double t280;
  double t281;
  double t283;
  double t284;
  double t285;
  double t148;
  double t150;
  double t155;
  double t156;
  double t157;
  double t302;
  double t303;
  double t229;
  double t236;
  double t305;
  double t307;
  double t308;
  double t309;
  double t310;
  double t311;
  double t312;
  double t241;
  double t243;
  double t244;
  double t246;
  double t247;
  double t248;
  double t249;
  double t250;
  double t253;
  double t256;
  double t257;
  double t258;
  double t261;
  double t55;
  double t62;
  double t73;
  double t76;
  double t77;
  double t87;
  double t295;
  double t296;
  double t298;
  double t330;
  double t331;
  double t89;
  double t183;
  double t186;
  double t191;
  double t201;
  double t203;
  double t357;
  double t358;
  double t359;
  double t360;
  double t361;
  double t352;
  double t382;
  double t383;
  double t291;
  double t167;
  double t362;
  double t363;
  double t398;
  double t399;
  double t366;
  double t367;
  double t368;
  double t369;
  double t370;
  double t371;
  double t372;
  double t373;
  double t390;
  double t391;
  double t282;
  double t286;
  double t415;
  double t416;
  double t417;
  double t418;
  double t419;
  double t289;
  double t290;
  double t292;
  double t293;
  double t294;
  double t299;
  double t136;
  double t137;
  double t138;
  double t139;
  double t153;
  double t158;
  double t161;
  double t164;
  double t166;
  double t169;
  double t170;
  double t171;
  double t172;
  double t175;
  double t384;
  double t385;
  double t436;
  double t437;
  double t388;
  double t389;
  double t392;
  double t393;
  double t394;
  double t395;
  double t332;
  double t333;
  double t334;
  double t335;
  double t336;
  double t337;
  double t338;
  double t339;
  double t340;
  double t341;
  double t342;
  double t343;
  double t346;
  double t407;
  double t408;
  double t409;
  double t410;
  double t411;
  double t412;
  double t413;
  double t64;
  double t78;
  double t92;
  double t93;
  double t101;
  double t113;
  double t115;
  double t430;
  double t431;
  double t320;
  double t321;
  double t322;
  double t323;
  double t324;
  double t325;
  double t326;
  double t428;
  double t429;
  double t432;
  double t433;
  double t434;
  double t448;
  double t449;
  double t450;
  double t194;
  double t204;
  double t205;
  double t207;
  double t208;
  double t210;
  double t211;
  double t213;
  double t214;
  double t215;
  double t237;
  double t238;
  double t262;
  double t263;
  double t265;
  double t468;
  double t469;
  double t470;
  double t471;
  double t348;
  double t349;
  double t350;
  double t351;
  double t353;
  double t354;
  double t355;
  double t356;
  double t364;
  double t365;
  double t374;
  double t375;
  double t376;
  double t480;
  double t490;
  double t491;
  double t492;
  t34 = Cos(var1[3]);
  t43 = Sin(var1[2]);
  t16 = Cos(var1[2]);
  t44 = Sin(var1[3]);
  t53 = Cos(var1[4]);
  t56 = Sin(var1[4]);
  t42 = t16*t34;
  t46 = t43*t44;
  t48 = t42 + t46;
  t65 = -1.*t34*t43;
  t66 = t16*t44;
  t70 = t65 + t66;
  t119 = 0.0265*t34;
  t120 = -0.0695*t44;
  t121 = t119 + t120;
  t123 = t16*t121;
  t124 = -0.0695*t34;
  t125 = -0.0265*t44;
  t126 = t124 + t125;
  t132 = -1.*t43*t126;
  t60 = -0.2375*t56;
  t133 = t34*t43;
  t134 = -1.*t16*t44;
  t135 = t133 + t134;
  t140 = -1.*t53;
  t141 = 1. + t140;
  t82 = t53*t48;
  t100 = -1.*t48*t56;
  t177 = Cos(var1[5]);
  t179 = Sin(var1[5]);
  t182 = Cos(var1[6]);
  t184 = Sin(var1[6]);
  t178 = t16*t177;
  t180 = -1.*t43*t179;
  t181 = t178 + t180;
  t195 = -1.*t177*t43;
  t196 = -1.*t16*t179;
  t199 = t195 + t196;
  t216 = -0.0265*t177;
  t217 = -0.0695*t179;
  t218 = t216 + t217;
  t219 = t16*t218;
  t220 = 0.0695*t177;
  t221 = -0.0265*t179;
  t225 = t220 + t221;
  t228 = -1.*t43*t225;
  t202 = -0.0265*t184;
  t231 = -1.*t16*t177;
  t233 = t43*t179;
  t235 = t231 + t233;
  t239 = -1.*t182;
  t240 = 1. + t239;
  t206 = -1.*t199*t184;
  t209 = t182*t199;
  t271 = -1.*t34;
  t272 = 1. + t271;
  t274 = -0.0265*t272;
  t275 = t274 + t120;
  t276 = -1.*t16*t275;
  t277 = -0.0695*t272;
  t278 = 0.0265*t44;
  t280 = t277 + t278;
  t281 = -1.*t43*t280;
  t283 = -1.*t16*t34;
  t284 = -1.*t43*t44;
  t285 = t283 + t284;
  t148 = -0.0265*t141;
  t150 = t148 + t60;
  t155 = -0.2375*t141;
  t156 = 0.0265*t56;
  t157 = t155 + t156;
  t302 = -1.*t177;
  t303 = 1. + t302;
  t229 = -0.15232*t199;
  t236 = -0.025413*t235;
  t305 = -0.0695*t303;
  t307 = t305 + t221;
  t308 = -1.*t43*t307;
  t309 = -0.0265*t303;
  t310 = 0.0695*t179;
  t311 = t309 + t310;
  t312 = -1.*t16*t311;
  t241 = -0.2375*t240;
  t243 = t241 + t202;
  t244 = t199*t243;
  t246 = -0.0265*t240;
  t247 = 0.2375*t184;
  t248 = t246 + t247;
  t249 = t235*t248;
  t250 = t182*t235;
  t253 = t250 + t206;
  t256 = -0.025226*t253;
  t257 = t235*t184;
  t258 = t209 + t257;
  t261 = -0.314506*t258;
  t55 = 0.0265*t53;
  t62 = t55 + t60;
  t73 = -0.2375*t53;
  t76 = -0.0265*t56;
  t77 = t73 + t76;
  t87 = -1.*t70*t56;
  t295 = t53*t70;
  t296 = -1.*t285*t56;
  t298 = t295 + t296;
  t330 = -1.*t43*t121;
  t331 = -1.*t16*t126;
  t89 = t82 + t87;
  t183 = -0.0265*t182;
  t186 = -0.2375*t184;
  t191 = t183 + t186;
  t201 = 0.2375*t182;
  t203 = t201 + t202;
  t357 = -1.*t43*t218;
  t358 = -1.*t16*t225;
  t359 = t177*t43;
  t360 = t16*t179;
  t361 = t359 + t360;
  t352 = -1.*t235*t184;
  t382 = t43*t275;
  t383 = -1.*t16*t280;
  t291 = t53*t285;
  t167 = t53*t135;
  t362 = -0.025413*t361;
  t363 = -0.15232*t235;
  t398 = -1.*t16*t307;
  t399 = t43*t311;
  t366 = t235*t243;
  t367 = t361*t248;
  t368 = t361*t184;
  t369 = t250 + t368;
  t370 = -0.314506*t369;
  t371 = t182*t361;
  t372 = t371 + t352;
  t373 = -0.025226*t372;
  t390 = -1.*t135*t56;
  t391 = t291 + t390;
  t282 = -0.15232*t70;
  t286 = -0.025367*t285;
  t415 = t43*t126;
  t416 = -0.0265*t34;
  t417 = 0.0695*t44;
  t418 = t416 + t417;
  t419 = t16*t418;
  t289 = t285*t150;
  t290 = t70*t157;
  t292 = t70*t56;
  t293 = t291 + t292;
  t294 = -0.025229*t293;
  t299 = -0.314514*t298;
  t136 = -0.15232*t135;
  t137 = -0.025367*t48;
  t138 = t123 + t132 + t136 + t137;
  t139 = 0.69051*t138;
  t153 = t48*t150;
  t158 = t135*t157;
  t161 = t135*t56;
  t164 = t82 + t161;
  t166 = -0.025229*t164;
  t169 = t167 + t100;
  t170 = -0.314514*t169;
  t171 = t123 + t132 + t153 + t158 + t166 + t170;
  t172 = 0.19964*t171;
  t175 = t139 + t172;
  t384 = -0.025367*t135;
  t385 = -0.15232*t285;
  t436 = t16*t126;
  t437 = -1.*t43*t418;
  t388 = t135*t150;
  t389 = t285*t157;
  t392 = -0.314514*t391;
  t393 = t285*t56;
  t394 = t167 + t393;
  t395 = -0.025229*t394;
  t332 = -0.025367*t70;
  t333 = -0.15232*t48;
  t334 = t330 + t331 + t332 + t333;
  t335 = 0.69051*t334;
  t336 = t70*t150;
  t337 = t48*t157;
  t338 = -0.314514*t89;
  t339 = t48*t56;
  t340 = t295 + t339;
  t341 = -0.025229*t340;
  t342 = t330 + t331 + t336 + t337 + t338 + t341;
  t343 = 0.19964*t342;
  t346 = t335 + t343;
  t407 = t285*t62;
  t408 = t135*t77;
  t409 = -0.025229*t391;
  t410 = -1.*t53*t135;
  t411 = t410 + t296;
  t412 = -0.314514*t411;
  t413 = t407 + t408 + t409 + t412;
  t64 = t48*t62;
  t78 = t70*t77;
  t92 = -0.025229*t89;
  t93 = -1.*t53*t70;
  t101 = t93 + t100;
  t113 = -0.314514*t101;
  t115 = t64 + t78 + t92 + t113;
  t430 = -1.*t53*t48;
  t431 = t430 + t390;
  t320 = t70*t62;
  t321 = t285*t77;
  t322 = -1.*t53*t285;
  t323 = t322 + t87;
  t324 = -0.314514*t323;
  t325 = -0.025229*t298;
  t326 = t320 + t321 + t324 + t325;
  t428 = t135*t62;
  t429 = t48*t77;
  t432 = -0.314514*t431;
  t433 = -0.025229*t169;
  t434 = t428 + t429 + t432 + t433;
  t448 = -0.0265*t53;
  t449 = 0.2375*t56;
  t450 = t448 + t449;
  t194 = t181*t191;
  t204 = t199*t203;
  t205 = -1.*t182*t181;
  t207 = t205 + t206;
  t208 = -0.025226*t207;
  t210 = -1.*t181*t184;
  t211 = t209 + t210;
  t213 = -0.314506*t211;
  t214 = t194 + t204 + t208 + t213;
  t215 = 0.05982798405705895*var2[6]*t214;
  t237 = t219 + t228 + t229 + t236;
  t238 = 0.69051*t237;
  t262 = t219 + t228 + t244 + t249 + t256 + t261;
  t263 = 0.19964*t262;
  t265 = t238 + t263;
  t468 = -0.0695*t177;
  t469 = 0.0265*t179;
  t470 = t468 + t469;
  t471 = t43*t470;
  t348 = t199*t191;
  t349 = t235*t203;
  t350 = -0.314506*t253;
  t351 = -1.*t182*t199;
  t353 = t351 + t352;
  t354 = -0.025226*t353;
  t355 = t348 + t349 + t350 + t354;
  t356 = 0.05982798405705895*var2[6]*t355;
  t364 = t357 + t358 + t362 + t363;
  t365 = 0.69051*t364;
  t374 = t357 + t358 + t366 + t367 + t370 + t373;
  t375 = 0.19964*t374;
  t376 = t365 + t375;
  t480 = t16*t470;
  t490 = -0.2375*t182;
  t491 = 0.0265*t184;
  t492 = t490 + t491;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=t215 + 0.2996793431028799*(0.69051*(t276 + t281 + t282 + t286) + 0.19964*(t276 + t281 + t289 + t290 + t294 + t299) + 0.69051*(t229 + t236 + t308 + t312) + 0.19964*(t244 + t249 + t256 + t261 + t308 + t312) + 1.5566*(0.026461*t16 + 0.046589*t43))*var2[2] + 0.2996793431028799*t175*var2[3] + 0.05982798405705895*t115*var2[4] + 0.2996793431028799*t265*var2[5];
  p_output1[7]=0;
  p_output1[8]=t356 + 0.2996793431028799*(0.69051*(t382 + t383 + t384 + t385) + 0.19964*(t382 + t383 + t388 + t389 + t392 + t395) + 0.69051*(t362 + t363 + t398 + t399) + 0.19964*(t366 + t367 + t370 + t373 + t398 + t399) + 1.5566*(0.046589*t16 - 0.026461*t43))*var2[2] + 0.2996793431028799*t346*var2[3] + 0.05982798405705895*t326*var2[4] + 0.2996793431028799*t376*var2[5];
  p_output1[9]=0.2996793431028799*t175*var2[2] + 0.2996793431028799*(0.69051*(t282 + t286 + t415 + t419) + 0.19964*(t289 + t290 + t294 + t299 + t415 + t419))*var2[3] + 0.05982798405705895*t413*var2[4];
  p_output1[10]=0;
  p_output1[11]=0.2996793431028799*t346*var2[2] + 0.2996793431028799*(0.69051*(t384 + t385 + t436 + t437) + 0.19964*(t388 + t389 + t392 + t395 + t436 + t437))*var2[3] + 0.05982798405705895*t434*var2[4];
  p_output1[12]=0.05982798405705895*t115*var2[2] + 0.05982798405705895*t413*var2[3] + 0.05982798405705895*(t408 - 0.314514*(t339 + t410) - 0.025229*t431 + t450*t48)*var2[4];
  p_output1[13]=0;
  p_output1[14]=0.05982798405705895*t326*var2[2] + 0.05982798405705895*t434*var2[3] + 0.05982798405705895*(-0.025229*t101 + t429 - 0.314514*(t292 + t430) + t450*t70)*var2[4];
  p_output1[15]=t215 + 0.2996793431028799*t265*var2[2] + 0.2996793431028799*(0.69051*(t219 + t229 + t236 + t471) + 0.19964*(t219 + t244 + t249 + t256 + t261 + t471))*var2[5];
  p_output1[16]=0;
  p_output1[17]=t356 + 0.2996793431028799*t376*var2[2] + 0.2996793431028799*(0.69051*(t357 + t362 + t363 + t480) + 0.19964*(t357 + t366 + t367 + t370 + t373 + t480))*var2[5];
  p_output1[18]=0.05982798405705895*t214*var2[2] + 0.05982798405705895*t214*var2[5] + 0.05982798405705895*(t194 - 0.314506*(t210 - 1.*t182*t361) - 0.025226*(t205 + t368) + t361*t492)*var2[6];
  p_output1[19]=0;
  p_output1[20]=0.05982798405705895*t355*var2[2] + 0.05982798405705895*t355*var2[5] + 0.05982798405705895*(-0.314506*t207 + t348 - 0.025226*(t181*t184 + t351) + t181*t492)*var2[6];
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 3, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "dJ_COMPosition.hh"

namespace SymFunction
{

void dJ_COMPosition_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
