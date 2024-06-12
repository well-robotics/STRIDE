/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:28 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
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


#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t26;
  double t14;
  double t15;
  double t27;
  double t34;
  double t10;
  double t35;
  double t37;
  double t38;
  double t16;
  double t30;
  double t32;
  double t33;
  double t40;
  double t41;
  double t100;
  double t112;
  double t113;
  double t43;
  double t45;
  double t47;
  double t48;
  double t49;
  double t59;
  double t60;
  double t63;
  double t65;
  double t66;
  double t67;
  double t68;
  double t70;
  double t71;
  double t75;
  double t81;
  double t93;
  double t117;
  double t123;
  double t124;
  double t137;
  double t138;
  double t171;
  double t177;
  double t180;
  double t181;
  double t182;
  double t188;
  double t197;
  double t189;
  double t201;
  double t206;
  double t209;
  double t210;
  double t211;
  double t212;
  double t218;
  double t219;
  double t220;
  double t208;
  double t213;
  double t235;
  double t236;
  double t237;
  double t214;
  double t175;
  double t178;
  double t179;
  double t183;
  double t184;
  double t185;
  double t186;
  double t187;
  double t243;
  double t244;
  double t247;
  double t248;
  double t249;
  double t250;
  double t251;
  double t159;
  double t163;
  double t149;
  double t152;
  double t258;
  double t259;
  double t260;
  double t261;
  double t262;
  double t263;
  double t264;
  double t267;
  double t268;
  double t269;
  double t270;
  double t271;
  double t272;
  double t273;
  double t274;
  double t275;
  double t233;
  double t234;
  double t215;
  double t216;
  double t297;
  double t288;
  double t289;
  double t290;
  double t291;
  double t292;
  double t293;
  double t294;
  double t295;
  double t296;
  double t310;
  double t311;
  double t312;
  double t313;
  double t314;
  double t315;
  double t328;
  double t329;
  double t330;
  double t331;
  double t332;
  double t333;
  double t334;
  double t335;
  double t337;
  double t338;
  double t339;
  double t343;
  double t344;
  double t345;
  double t336;
  double t340;
  double t341;
  double t342;
  double t347;
  double t348;
  double t352;
  double t353;
  double t354;
  double t355;
  double t364;
  double t365;
  double t357;
  double t367;
  double t368;
  double t359;
  double t322;
  double t323;
  double t324;
  double t325;
  double t326;
  double t327;
  double t387;
  double t388;
  double t389;
  double t390;
  double t391;
  double t392;
  double t393;
  double t394;
  double t396;
  double t397;
  double t398;
  double t381;
  double t382;
  double t383;
  double t384;
  double t385;
  double t386;
  double t395;
  double t399;
  double t400;
  double t402;
  double t403;
  double t404;
  double t409;
  double t410;
  double t411;
  double t408;
  double t413;
  double t414;
  double t418;
  double t427;
  double t428;
  double t420;
  double t430;
  double t431;
  double t422;
  double t443;
  double t444;
  double t445;
  double t446;
  double t448;
  double t449;
  double t450;
  double t451;
  double t455;
  double t456;
  double t476;
  double t477;
  double t478;
  double t479;
  double t481;
  double t482;
  double t483;
  double t484;
  double t488;
  double t489;
  t26 = Cos(var1[3]);
  t14 = Cos(var1[4]);
  t15 = Sin(var1[3]);
  t27 = Sin(var1[4]);
  t34 = Cos(var1[2]);
  t10 = Sin(var1[2]);
  t35 = t26*t14;
  t37 = -1.*t15*t27;
  t38 = t35 + t37;
  t16 = -1.*t14*t15;
  t30 = -1.*t26*t27;
  t32 = t16 + t30;
  t33 = t10*t32;
  t40 = t34*t38;
  t41 = t33 + t40;
  t100 = t34*t26;
  t112 = -1.*t10*t15;
  t113 = t100 + t112;
  t43 = t14*t15;
  t45 = t26*t27;
  t47 = t43 + t45;
  t48 = t34*t47;
  t49 = t10*t38;
  t59 = t48 + t49;
  t60 = 1.70432*t41*t59;
  t63 = t34*t32;
  t65 = -1.*t26*t14;
  t66 = t15*t27;
  t67 = t65 + t66;
  t68 = t10*t67;
  t70 = t63 + t68;
  t71 = 1.70432*t41*t70;
  t75 = -1.*t26*t10;
  t81 = -1.*t34*t15;
  t93 = t75 + t81;
  t117 = 6.8522*t93*t113;
  t123 = t26*t10;
  t124 = t34*t15;
  t137 = t123 + t124;
  t138 = 6.8522*t137*t113;
  t171 = Cos(var1[5]);
  t177 = Sin(var1[5]);
  t180 = t34*t171;
  t181 = -1.*t10*t177;
  t182 = t180 + t181;
  t188 = Cos(var1[6]);
  t197 = Sin(var1[6]);
  t189 = -1.*t188*t177;
  t201 = -1.*t171*t197;
  t206 = t189 + t201;
  t209 = t171*t188;
  t210 = -1.*t177*t197;
  t211 = t209 + t210;
  t212 = t34*t211;
  t218 = t188*t177;
  t219 = t171*t197;
  t220 = t218 + t219;
  t208 = t10*t206;
  t213 = t208 + t212;
  t235 = t34*t220;
  t236 = t10*t211;
  t237 = t235 + t236;
  t214 = t34*t206;
  t175 = -1.*t171*t10;
  t178 = -1.*t34*t177;
  t179 = t175 + t178;
  t183 = 6.8522*t179*t182;
  t184 = t171*t10;
  t185 = t34*t177;
  t186 = t184 + t185;
  t187 = 6.8522*t186*t182;
  t243 = 1.70432*t213*t237;
  t244 = -1.*t171*t188;
  t247 = t177*t197;
  t248 = t244 + t247;
  t249 = t10*t248;
  t250 = t214 + t249;
  t251 = 1.70432*t213*t250;
  t159 = -1.*t10*t47;
  t163 = t159 + t40;
  t149 = -1.*t10*t38;
  t152 = t63 + t149;
  t258 = 0.85216*t41*t163;
  t259 = 0.85216*t152*t59;
  t260 = -1.*t10*t32;
  t261 = t34*t67;
  t262 = t260 + t261;
  t263 = 0.85216*t41*t262;
  t264 = 0.85216*t152*t70;
  t267 = Power(t93,2);
  t268 = 3.4261*t267;
  t269 = 3.4261*t93*t137;
  t270 = Power(t113,2);
  t271 = 3.4261*t270;
  t272 = -1.*t34*t26;
  t273 = t10*t15;
  t274 = t272 + t273;
  t275 = 3.4261*t113*t274;
  t233 = -1.*t10*t220;
  t234 = t233 + t212;
  t215 = -1.*t10*t211;
  t216 = t214 + t215;
  t297 = -1.*t10*t206;
  t288 = Power(t179,2);
  t289 = 3.4261*t288;
  t290 = 3.4261*t179*t186;
  t291 = Power(t182,2);
  t292 = 3.4261*t291;
  t293 = -1.*t34*t171;
  t294 = t10*t177;
  t295 = t293 + t294;
  t296 = 3.4261*t182*t295;
  t310 = 0.85216*t213*t234;
  t311 = 0.85216*t216*t237;
  t312 = t34*t248;
  t313 = t297 + t312;
  t314 = 0.85216*t213*t313;
  t315 = 0.85216*t216*t250;
  t328 = -1.*t14;
  t329 = 1. + t328;
  t330 = 0.5*t329;
  t331 = 0.671885*t14;
  t332 = t330 + t331;
  t333 = t332*t15;
  t334 = 0.171885*t26*t27;
  t335 = t333 + t334;
  t337 = t26*t332;
  t338 = -0.171885*t15*t27;
  t339 = t337 + t338;
  t343 = -1.*t332*t15;
  t344 = -0.171885*t26*t27;
  t345 = t343 + t344;
  t336 = -1.*t335*t38;
  t340 = -1.*t32*t339;
  t341 = t336 + t340;
  t342 = 0.85216*t41*t341;
  t347 = t335*t38;
  t348 = t32*t339;
  t352 = t335*t47;
  t353 = t38*t339;
  t354 = t352 + t353;
  t355 = 0.85216*t354*t70;
  t364 = -0.171885*t14*t15;
  t365 = t364 + t344;
  t357 = -1.*t32*t335;
  t367 = 0.171885*t26*t14;
  t368 = t367 + t338;
  t359 = -1.*t339*t67;
  t322 = Power(t26,2);
  t323 = 0.1494*t322;
  t324 = Power(t15,2);
  t325 = 0.1494*t324;
  t326 = t323 + t325;
  t327 = 3.4261*t93*t326;
  t387 = -1.*t188;
  t388 = 1. + t387;
  t389 = 0.5*t388;
  t390 = 0.671885*t188;
  t391 = t389 + t390;
  t392 = t391*t177;
  t393 = 0.171885*t171*t197;
  t394 = t392 + t393;
  t396 = t171*t391;
  t397 = -0.171885*t177*t197;
  t398 = t396 + t397;
  t381 = Power(t171,2);
  t382 = 0.1494*t381;
  t383 = Power(t177,2);
  t384 = 0.1494*t383;
  t385 = t382 + t384;
  t386 = 3.4261*t179*t385;
  t395 = -1.*t394*t211;
  t399 = -1.*t206*t398;
  t400 = t395 + t399;
  t402 = t394*t220;
  t403 = t211*t398;
  t404 = t402 + t403;
  t409 = -1.*t391*t177;
  t410 = -0.171885*t171*t197;
  t411 = t409 + t410;
  t408 = 0.85216*t213*t400;
  t413 = t394*t211;
  t414 = t206*t398;
  t418 = 0.85216*t404*t250;
  t427 = -0.171885*t188*t177;
  t428 = t427 + t410;
  t420 = -1.*t206*t394;
  t430 = 0.171885*t171*t188;
  t431 = t430 + t397;
  t422 = -1.*t398*t248;
  t443 = 0.51185934*t93;
  t444 = t332*t27;
  t445 = -0.171885*t14*t27;
  t446 = t444 + t445;
  t448 = t332*t14;
  t449 = Power(t27,2);
  t450 = 0.171885*t449;
  t451 = t448 + t450;
  t455 = 0.85216*t446*t41;
  t456 = 0.85216*t451*t70;
  t476 = 0.51185934*t179;
  t477 = t391*t197;
  t478 = -0.171885*t188*t197;
  t479 = t477 + t478;
  t481 = t391*t188;
  t482 = Power(t197,2);
  t483 = 0.171885*t482;
  t484 = t481 + t483;
  t488 = 0.85216*t479*t213;
  t489 = 0.85216*t484*t250;
  p_output1[0]=var2[0]*(-0.5*(t117 + t138 + t183 + t187 + 1.70432*t213*t216 + 1.70432*t234*t237 + 1.70432*t152*t41 + 1.70432*t163*t59)*var2[2] - 0.5*(t117 + t138 + t60 + t71)*var2[3] - 0.5*(t60 + t71)*var2[4] - 0.5*(t183 + t187 + t243 + t251)*var2[5] - 0.5*(t243 + t251)*var2[6]);
  p_output1[1]=var2[0]*(-0.5*(0.85216*Power(t152,2) + 0.85216*Power(t163,2) + 0.85216*Power(t216,2) + 0.85216*Power(t234,2) + t268 + t269 + t271 + t275 + t289 + t290 + t292 + t296 + 0.85216*t213*(t297 - 1.*t211*t34) + 0.85216*t237*(t215 - 1.*t220*t34) + 0.85216*(t260 - 1.*t34*t38)*t41 + 0.85216*(t149 - 1.*t34*t47)*t59)*var2[2] - 0.5*(t258 + t259 + t263 + t264 + t268 + t269 + t271 + t275)*var2[3] - 0.5*(t258 + t259 + t263 + t264)*var2[4] - 0.5*(t289 + t290 + t292 + t296 + t310 + t311 + t314 + t315)*var2[5] - 0.5*(t310 + t311 + t314 + t315)*var2[6]);
  p_output1[2]=var2[0]*(-0.5*(-3.70591*t10 + t327 + 0.85216*t163*t341 + 0.85216*t152*t354 + t386 + 0.85216*t234*t400 + 0.85216*t216*t404)*var2[2] - 0.5*(t327 + t342 + t355 + 0.85216*t41*(t347 + t348 + t345*t38 + t339*t47) + 0.85216*(-1.*t32*t345 + t357 + t359 - 1.*t339*t38)*t59)*var2[3] - 0.5*(t342 + t355 + 0.85216*t41*(t347 + t348 + t365*t38 + t368*t47) + 0.85216*(t357 + t359 - 1.*t32*t365 - 1.*t368*t38)*t59)*var2[4] - 0.5*(t386 + t408 + 0.85216*t213*(t220*t398 + t211*t411 + t413 + t414) + t418 + 0.85216*t237*(-1.*t211*t398 - 1.*t206*t411 + t420 + t422))*var2[5] - 0.5*(t408 + t418 + 0.85216*t237*(t420 + t422 - 1.*t206*t428 - 1.*t211*t431) + 0.85216*t213*(t413 + t414 + t211*t428 + t220*t431))*var2[6]);
  p_output1[3]=var2[0]*(-0.5*(t443 + 0.85216*t163*t446 + 0.85216*t152*t451)*var2[2] - 0.5*(t443 + t455 + t456)*var2[3] - 0.5*(0.85216*(0.171885*t14*t27 - 1.*t27*t332)*t41 + t455 + t456 + 0.85216*(-0.171885*Power(t14,2) + t448)*t59)*var2[4]);
  p_output1[4]=var2[0]*(-0.0732367608*t152*var2[2] - 0.0732367608*t70*var2[3] - 0.0732367608*t70*var2[4]);
  p_output1[5]=var2[0]*(-0.5*(t476 + 0.85216*t234*t479 + 0.85216*t216*t484)*var2[2] - 0.5*(t476 + t488 + t489)*var2[5] - 0.5*(0.85216*t213*(0.171885*t188*t197 - 1.*t197*t391) + 0.85216*t237*(-0.171885*Power(t188,2) + t481) + t488 + t489)*var2[6]);
  p_output1[6]=var2[0]*(-0.0732367608*t216*var2[2] - 0.0732367608*t250*var2[5] - 0.0732367608*t250*var2[6]);
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 1, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1,var2);


}

#else // MATLAB_MEX_FILE

#include "Ce1_vec1_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
