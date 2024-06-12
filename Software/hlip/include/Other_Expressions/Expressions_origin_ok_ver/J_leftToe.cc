/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:28 GMT-05:00
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
  double t92;
  double t113;
  double t125;
  double t136;
  double t78;
  double t161;
  double t194;
  double t202;
  double t203;
  double t207;
  double t175;
  double t186;
  double t191;
  double t213;
  double t214;
  double t215;
  double t133;
  double t137;
  double t148;
  double t166;
  double t167;
  double t169;
  double t204;
  double t208;
  double t209;
  double t217;
  double t220;
  double t221;
  double t253;
  double t256;
  double t257;
  double t241;
  double t211;
  double t229;
  double t233;
  double t236;
  double t237;
  double t238;
  double t243;
  double t244;
  double t246;
  double t270;
  double t274;
  double t277;
  double t282;
  double t286;
  double t250;
  double t258;
  double t261;
  double t262;
  double t263;
  double t265;
  double t266;
  double t267;
  double t268;
  double t292;
  double t293;
  double t294;
  double t295;
  double t298;
  double t299;
  double t301;
  double t305;
  t92 = Cos(var1[5]);
  t113 = -1.*t92;
  t125 = 1. + t113;
  t136 = Sin(var1[5]);
  t78 = Cos(var1[2]);
  t161 = Sin(var1[2]);
  t194 = Cos(var1[6]);
  t202 = -1.*t194;
  t203 = 1. + t202;
  t207 = Sin(var1[6]);
  t175 = t78*t92;
  t186 = -1.*t161*t136;
  t191 = t175 + t186;
  t213 = -1.*t92*t161;
  t214 = -1.*t78*t136;
  t215 = t213 + t214;
  t133 = -0.0695*t125;
  t137 = -0.0265*t136;
  t148 = t133 + t137;
  t166 = -0.0265*t125;
  t167 = 0.0695*t136;
  t169 = t166 + t167;
  t204 = -0.2375*t203;
  t208 = -0.0265*t207;
  t209 = t204 + t208;
  t217 = -0.0265*t203;
  t220 = 0.2375*t207;
  t221 = t217 + t220;
  t253 = -1.*t78*t92;
  t256 = t161*t136;
  t257 = t253 + t256;
  t241 = t194*t215;
  t211 = t191*t209;
  t229 = t215*t221;
  t233 = t194*t191;
  t236 = t215*t207;
  t237 = t233 + t236;
  t238 = 0.0225*t237;
  t243 = -1.*t191*t207;
  t244 = t241 + t243;
  t246 = -0.0265*t244;
  t270 = -0.0265*t92;
  t274 = -0.0695*t136;
  t277 = t270 + t274;
  t282 = 0.0695*t92;
  t286 = t282 + t137;
  t250 = t215*t209;
  t258 = t257*t221;
  t261 = t194*t257;
  t262 = -1.*t215*t207;
  t263 = t261 + t262;
  t265 = -0.0265*t263;
  t266 = t257*t207;
  t267 = t241 + t266;
  t268 = 0.0225*t267;
  t292 = t92*t161;
  t293 = t78*t136;
  t294 = t292 + t293;
  t295 = -0.0265*t194;
  t298 = -0.2375*t207;
  t299 = t295 + t298;
  t301 = 0.2375*t194;
  t305 = t301 + t208;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=1.;
  p_output1[4]=-1.*t161*t169 + t211 + t229 + t238 + t246 + t148*t78;
  p_output1[5]=-1.*t148*t161 + t250 + t258 + t265 + t268 - 1.*t169*t78;
  p_output1[6]=0;
  p_output1[7]=0;
  p_output1[8]=0;
  p_output1[9]=0;
  p_output1[10]=t211 + t229 + t238 + t246 + t161*t277 + t286*t78;
  p_output1[11]=t250 + t258 + t265 + t268 - 1.*t161*t286 + t277*t78;
  p_output1[12]=-0.0265*(t243 - 1.*t194*t294) + 0.0225*(t233 - 1.*t207*t294) + t294*t299 + t191*t305;
  p_output1[13]=0.0225*t244 - 0.0265*(-1.*t191*t194 + t262) + t191*t299 + t215*t305;
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
  plhs[0] = mxCreateDoubleMatrix((mwSize) 2, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "J_leftToe.hh"

namespace SymFunction
{

void J_leftToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
