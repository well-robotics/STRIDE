/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:54 GMT-05:00
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
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t199;
  double t169;
  double t176;
  double t203;
  double t158;
  double t222;
  double t223;
  double t224;
  double t225;
  double t232;
  double t235;
  double t238;
  double t226;
  double t191;
  double t204;
  double t209;
  double t249;
  double t210;
  double t227;
  double t244;
  double t246;
  double t248;
  double t250;
  double t255;
  double t256;
  double t257;
  double t265;
  double t266;
  double t268;
  double t261;
  double t262;
  double t263;
  double t269;
  double t270;
  double t271;
  double t274;
  double t275;
  double t278;
  double t279;
  double t253;
  double t254;
  double t239;
  double t241;
  double t243;
  double t293;
  double t296;
  double t297;
  double t308;
  double t309;
  double t302;
  double t303;
  double t313;
  double t315;
  double t316;
  double t317;
  double t319;
  double t320;
  double t321;
  double t322;
  double t326;
  double t327;
  double t328;
  double t330;
  double t331;
  double t305;
  double t307;
  double t318;
  double t323;
  double t324;
  double t339;
  double t340;
  double t341;
  double t325;
  double t333;
  double t335;
  double t349;
  double t354;
  double t355;
  double t356;
  double t358;
  double t346;
  double t347;
  double t363;
  double t345;
  double t360;
  double t375;
  double t376;
  double t377;
  double t379;
  double t380;
  double t366;
  double t368;
  double t384;
  double t385;
  double t386;
  double t378;
  double t381;
  double t382;
  double t370;
  double t399;
  double t400;
  double t401;
  double t402;
  double t403;
  double t405;
  double t406;
  double t407;
  double t408;
  double t409;
  double t413;
  double t414;
  double t437;
  double t438;
  double t439;
  t199 = Cos(var1[3]);
  t169 = Cos(var1[4]);
  t176 = Sin(var1[3]);
  t203 = Sin(var1[4]);
  t158 = Sin(var1[2]);
  t222 = Cos(var1[2]);
  t223 = t199*t169;
  t224 = -1.*t176*t203;
  t225 = t223 + t224;
  t232 = t169*t176;
  t235 = t199*t203;
  t238 = t232 + t235;
  t226 = t222*t225;
  t191 = -1.*t169*t176;
  t204 = -1.*t199*t203;
  t209 = t191 + t204;
  t249 = -1.*t158*t225;
  t210 = -1.*t158*t209;
  t227 = t210 + t226;
  t244 = t158*t238;
  t246 = t244 + t226;
  t248 = t222*t238;
  t250 = t248 + t249;
  t255 = t222*t209;
  t256 = t158*t225;
  t257 = t255 + t256;
  t265 = -1.*t199*t169;
  t266 = t176*t203;
  t268 = t265 + t266;
  t261 = 0.19964*t227*t246;
  t262 = 0.19964*t250*t257;
  t263 = t158*t209;
  t269 = t222*t268;
  t270 = t263 + t269;
  t271 = 0.19964*t227*t270;
  t274 = -1.*t158*t268;
  t275 = t255 + t274;
  t278 = 0.19964*t257*t275;
  t279 = t261 + t262 + t271 + t278;
  t253 = -1.*t222*t209;
  t254 = t253 + t249;
  t239 = -1.*t158*t238;
  t241 = -1.*t222*t225;
  t243 = t239 + t241;
  t293 = 0.39928*t227*t250;
  t296 = 0.39928*t227*t275;
  t297 = t293 + t296;
  t308 = -1.*t169;
  t309 = 1. + t308;
  t302 = -1.*t199;
  t303 = 1. + t302;
  t313 = -0.2375*t309;
  t315 = -0.314514*t169;
  t316 = 0.0012709999999999978*t203;
  t317 = t313 + t315 + t316;
  t319 = -0.0265*t309;
  t320 = -0.025229*t169;
  t321 = 0.07701400000000003*t203;
  t322 = t319 + t320 + t321;
  t326 = -0.0265*t303;
  t327 = -0.0695*t176;
  t328 = -1.*t176*t317;
  t330 = t199*t322;
  t331 = t326 + t327 + t328 + t330;
  t305 = -0.0695*t303;
  t307 = 0.0265*t176;
  t318 = t199*t317;
  t323 = t176*t322;
  t324 = t305 + t307 + t318 + t323;
  t339 = t331*t209;
  t340 = t324*t225;
  t341 = t339 + t340;
  t325 = -1.*t324*t238;
  t333 = -1.*t331*t225;
  t335 = t325 + t333;
  t349 = -0.0695*t199;
  t354 = -0.0265*t176;
  t355 = -1.*t199*t317;
  t356 = -1.*t176*t322;
  t358 = t349 + t354 + t355 + t356;
  t346 = 0.0265*t199;
  t347 = t346 + t327 + t328 + t330;
  t363 = 0.19964*t227*t341;
  t345 = -1.*t331*t209;
  t360 = -1.*t324*t225;
  t375 = 0.07701400000000003*t169;
  t376 = -0.0012709999999999978*t203;
  t377 = t375 + t376;
  t379 = 0.0012709999999999978*t169;
  t380 = t379 + t321;
  t366 = 0.19964*t335*t275;
  t368 = t324*t209;
  t384 = t199*t377;
  t385 = -1.*t176*t380;
  t386 = t384 + t385;
  t378 = t176*t377;
  t381 = t199*t380;
  t382 = t378 + t381;
  t370 = t331*t268;
  t399 = -0.0695*t169;
  t400 = -1.*t169*t317;
  t401 = 0.0265*t203;
  t402 = t322*t203;
  t403 = t399 + t400 + t401 + t402;
  t405 = 0.0265*t169;
  t406 = t169*t322;
  t407 = 0.0695*t203;
  t408 = t317*t203;
  t409 = t405 + t406 + t407 + t408;
  t413 = 0.19964*t403*t227;
  t414 = 0.19964*t409*t275;
  t437 = 0.015375074960000006*t227;
  t438 = 0.0002537424399999996*t275;
  t439 = t437 + t438;
  p_output1[0]=var2[1]*(-0.5*(0.19964*Power(t227,2) + 0.19964*t243*t246 + 0.19964*Power(t250,2) + 0.19964*t254*t257)*var2[2] - 0.5*t279*var2[3] - 0.5*t279*var2[4]);
  p_output1[1]=var2[1]*(-0.5*(0.39928*t243*t250 + 0.39928*t227*t254)*var2[2] - 0.5*t297*var2[3] - 0.5*t297*var2[4]);
  p_output1[2]=var2[1]*(-0.5*(0.19964*t254*t335 + 0.19964*t243*t341)*var2[2] - 0.5*(0.19964*t227*(t345 - 1.*t238*t347 - 1.*t225*t358 + t360) + t363 + t366 + 0.19964*t250*(t225*t347 + t209*t358 + t368 + t370))*var2[3] - 0.5*(t363 + t366 + 0.19964*t250*(t368 + t370 + t225*t382 + t209*t386) + 0.19964*t227*(t345 + t360 - 1.*t238*t382 - 1.*t225*t386))*var2[4]);
  p_output1[3]=var2[1]*(-0.5*(0.19964*t243*t403 + 0.19964*t254*t409)*var2[2] - 0.5*(t413 + t414)*var2[3] - 0.5*(0.19964*t227*(0.0695*t169 - 0.0265*t203 + t169*t317 - 1.*t203*t322 + t169*t377 + t203*t380) + 0.19964*t250*(t203*t377 - 1.*t169*t380 + t405 + t406 + t407 + t408) + t413 + t414)*var2[4]);
  p_output1[4]=var2[1]*(-0.5*(0.015375074960000006*t243 + 0.0002537424399999996*t254)*var2[2] - 0.5*t439*var2[3] - 0.5*t439*var2[4]);
  p_output1[5]=0;
  p_output1[6]=0;
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

#include "Ce1_vec_L3_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L3_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
