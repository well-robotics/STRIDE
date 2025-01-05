/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:57 GMT-05:00
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
  double t231;
  double t226;
  double t227;
  double t232;
  double t258;
  double t259;
  double t254;
  double t255;
  double t260;
  double t261;
  double t262;
  double t263;
  double t269;
  double t270;
  double t271;
  double t275;
  double t247;
  double t248;
  double t250;
  double t246;
  double t282;
  double t283;
  double t284;
  double t223;
  double t287;
  double t290;
  double t291;
  double t293;
  double t296;
  double t230;
  double t235;
  double t243;
  double t256;
  double t257;
  double t266;
  double t278;
  double t280;
  double t321;
  double t323;
  double t326;
  double t313;
  double t315;
  double t316;
  double t285;
  double t298;
  double t299;
  double t251;
  double t339;
  double t340;
  double t341;
  double t342;
  double t343;
  double t330;
  double t335;
  double t361;
  double t362;
  double t363;
  double t354;
  double t327;
  double t346;
  double t381;
  double t383;
  double t384;
  double t388;
  double t389;
  double t359;
  double t366;
  double t367;
  double t368;
  double t369;
  double t370;
  double t372;
  double t394;
  double t395;
  double t396;
  double t387;
  double t391;
  double t392;
  double t374;
  double t305;
  double t244;
  double t252;
  double t336;
  double t344;
  double t348;
  double t301;
  double t307;
  double t371;
  double t373;
  double t375;
  double t425;
  double t393;
  double t397;
  double t398;
  double t426;
  double t427;
  double t428;
  double t401;
  double t402;
  double t403;
  double t449;
  double t450;
  double t451;
  double t452;
  double t453;
  double t455;
  double t456;
  double t457;
  double t458;
  double t459;
  t231 = Cos(var1[3]);
  t226 = Cos(var1[4]);
  t227 = Sin(var1[3]);
  t232 = Sin(var1[4]);
  t258 = -1.*t226;
  t259 = 1. + t258;
  t254 = -1.*t231;
  t255 = 1. + t254;
  t260 = -0.2375*t259;
  t261 = -0.314514*t226;
  t262 = 0.0012709999999999978*t232;
  t263 = t260 + t261 + t262;
  t269 = -0.0265*t259;
  t270 = -0.025229*t226;
  t271 = 0.07701400000000003*t232;
  t275 = t269 + t270 + t271;
  t247 = t231*t226;
  t248 = -1.*t227*t232;
  t250 = t247 + t248;
  t246 = Cos(var1[2]);
  t282 = t226*t227;
  t283 = t231*t232;
  t284 = t282 + t283;
  t223 = Sin(var1[2]);
  t287 = -0.0265*t255;
  t290 = -0.0695*t227;
  t291 = -1.*t227*t263;
  t293 = t231*t275;
  t296 = t287 + t290 + t291 + t293;
  t230 = -1.*t226*t227;
  t235 = -1.*t231*t232;
  t243 = t230 + t235;
  t256 = -0.0695*t255;
  t257 = 0.0265*t227;
  t266 = t231*t263;
  t278 = t227*t275;
  t280 = t256 + t257 + t266 + t278;
  t321 = t246*t243;
  t323 = t223*t250;
  t326 = t321 + t323;
  t313 = t296*t243;
  t315 = t280*t250;
  t316 = t313 + t315;
  t285 = -1.*t280*t284;
  t298 = -1.*t296*t250;
  t299 = t285 + t298;
  t251 = t246*t250;
  t339 = -0.0695*t231;
  t340 = -0.0265*t227;
  t341 = -1.*t231*t263;
  t342 = -1.*t227*t275;
  t343 = t339 + t340 + t341 + t342;
  t330 = 0.0265*t231;
  t335 = t330 + t290 + t291 + t293;
  t361 = -1.*t231*t226;
  t362 = t227*t232;
  t363 = t361 + t362;
  t354 = 0.19964*t326*t316;
  t327 = -1.*t296*t243;
  t346 = -1.*t280*t250;
  t381 = 0.07701400000000003*t226;
  t383 = -0.0012709999999999978*t232;
  t384 = t381 + t383;
  t388 = 0.0012709999999999978*t226;
  t389 = t388 + t271;
  t359 = t223*t243;
  t366 = t246*t363;
  t367 = t359 + t366;
  t368 = 0.19964*t299*t367;
  t369 = t223*t284;
  t370 = t369 + t251;
  t372 = t280*t243;
  t394 = t231*t384;
  t395 = -1.*t227*t389;
  t396 = t394 + t395;
  t387 = t227*t384;
  t391 = t231*t389;
  t392 = t387 + t391;
  t374 = t296*t363;
  t305 = -1.*t223*t250;
  t244 = -1.*t223*t243;
  t252 = t244 + t251;
  t336 = -1.*t335*t284;
  t344 = -1.*t343*t250;
  t348 = t327 + t336 + t344 + t346;
  t301 = t246*t284;
  t307 = t301 + t305;
  t371 = t343*t243;
  t373 = t335*t250;
  t375 = t371 + t372 + t373 + t374;
  t425 = 0.19964*t252*t316;
  t393 = -1.*t392*t284;
  t397 = -1.*t396*t250;
  t398 = t327 + t393 + t346 + t397;
  t426 = -1.*t223*t363;
  t427 = t321 + t426;
  t428 = 0.19964*t299*t427;
  t401 = t396*t243;
  t402 = t392*t250;
  t403 = t372 + t401 + t402 + t374;
  t449 = 0.0265*t226;
  t450 = t226*t275;
  t451 = 0.0695*t232;
  t452 = t263*t232;
  t453 = t449 + t450 + t451 + t452;
  t455 = -0.0695*t226;
  t456 = -1.*t226*t263;
  t457 = 0.0265*t232;
  t458 = t275*t232;
  t459 = t455 + t456 + t457 + t458;
  p_output1[0]=var2[2]*(-0.5*(0.19964*t252*t299 + 0.19964*t307*t316)*var2[2] - 0.5*(0.19964*t326*t348 + t354 + t368 + 0.19964*t370*t375)*var2[3] - 0.5*(t354 + t368 + 0.19964*t326*t398 + 0.19964*t370*t403)*var2[4]);
  p_output1[1]=var2[2]*(-0.5*(0.19964*t299*(-1.*t243*t246 + t305) + 0.19964*(-1.*t246*t250 - 1.*t223*t284)*t316)*var2[2] - 0.5*(0.19964*t252*t348 + 0.19964*t307*t375 + t425 + t428)*var2[3] - 0.5*(0.19964*t252*t398 + 0.19964*t307*t403 + t425 + t428)*var2[4]);
  p_output1[2]=var2[2]*(-0.5*(0.39928*t299*t348 + 0.39928*t316*t375)*var2[3] - 0.5*(0.39928*t299*t398 + 0.39928*t316*t403)*var2[4]);
  p_output1[3]=var2[2]*(-0.5*(0.19964*t348*t453 + 0.19964*t375*t459)*var2[3] - 0.5*(0.19964*t299*(0.0695*t226 - 0.0265*t232 + t226*t263 - 1.*t232*t275 + t226*t384 + t232*t389) + 0.19964*t316*(t232*t384 - 1.*t226*t389 + t449 + t450 + t451 + t452) + 0.19964*t398*t453 + 0.19964*t403*t459)*var2[4]);
  p_output1[4]=var2[2]*(-0.5*(0.0002537424399999996*t348 + 0.015375074960000006*t375)*var2[3] - 0.5*(0.0002537424399999996*t398 + 0.015375074960000006*t403)*var2[4]);
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

#include "Ce1_vec_L3_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L3_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
