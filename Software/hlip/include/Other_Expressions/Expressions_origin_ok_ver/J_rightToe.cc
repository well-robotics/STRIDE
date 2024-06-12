/*
 * Automatically Generated from Mathematica.
 * Wed 20 Mar 2024 20:03:30 GMT-05:00
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
  double t137;
  double t150;
  double t166;
  double t170;
  double t133;
  double t209;
  double t247;
  double t248;
  double t249;
  double t257;
  double t267;
  double t269;
  double t270;
  double t237;
  double t241;
  double t244;
  double t167;
  double t175;
  double t204;
  double t217;
  double t220;
  double t221;
  double t256;
  double t261;
  double t263;
  double t274;
  double t278;
  double t282;
  double t305;
  double t309;
  double t310;
  double t292;
  double t315;
  double t335;
  double t336;
  double t337;
  double t324;
  double t325;
  double t327;
  double t332;
  double t333;
  double t288;
  double t343;
  double t340;
  double t362;
  double t363;
  double t366;
  double t367;
  double t369;
  double t370;
  double t371;
  double t289;
  double t290;
  t137 = Cos(var1[3]);
  t150 = -1.*t137;
  t166 = 1. + t150;
  t170 = Sin(var1[3]);
  t133 = Sin(var1[2]);
  t209 = Cos(var1[2]);
  t247 = Cos(var1[4]);
  t248 = -1.*t247;
  t249 = 1. + t248;
  t257 = Sin(var1[4]);
  t267 = t209*t137;
  t269 = t133*t170;
  t270 = t267 + t269;
  t237 = -1.*t137*t133;
  t241 = t209*t170;
  t244 = t237 + t241;
  t167 = -0.0265*t166;
  t175 = -0.0695*t170;
  t204 = t167 + t175;
  t217 = -0.0695*t166;
  t220 = 0.0265*t170;
  t221 = t217 + t220;
  t256 = -0.0265*t249;
  t261 = -0.2375*t257;
  t263 = t256 + t261;
  t274 = -0.2375*t249;
  t278 = 0.0265*t257;
  t282 = t274 + t278;
  t305 = -1.*t209*t137;
  t309 = -1.*t133*t170;
  t310 = t305 + t309;
  t292 = t247*t244;
  t315 = t247*t310;
  t335 = t137*t133;
  t336 = -1.*t209*t170;
  t337 = t335 + t336;
  t324 = 0.0265*t137;
  t325 = t324 + t175;
  t327 = -0.0695*t137;
  t332 = -0.0265*t170;
  t333 = t327 + t332;
  t288 = t247*t270;
  t343 = t247*t337;
  t340 = -1.*t337*t257;
  t362 = -1.*t270*t257;
  t363 = t343 + t362;
  t366 = 0.0265*t247;
  t367 = t366 + t261;
  t369 = -0.2375*t247;
  t370 = -0.0265*t257;
  t371 = t369 + t370;
  t289 = -1.*t244*t257;
  t290 = t288 + t289;
  p_output1[0]=1.;
  p_output1[1]=0;
  p_output1[2]=0;
  p_output1[3]=1.;
  p_output1[4]=-1.*t133*t204 + t209*t221 + t244*t263 + t270*t282 + 0.0225*t290 - 0.0265*(t257*t270 + t292);
  p_output1[5]=-1.*t204*t209 - 1.*t133*t221 + t244*t282 + t263*t310 + 0.0225*(t292 - 1.*t257*t310) - 0.0265*(t244*t257 + t315);
  p_output1[6]=t282*t310 + t133*t325 + t209*t333 + t263*t337 + 0.0225*(t315 + t340) - 0.0265*(t257*t310 + t343);
  p_output1[7]=t263*t270 + t209*t325 - 1.*t133*t333 + t282*t337 - 0.0265*(t288 + t257*t337) + 0.0225*t363;
  p_output1[8]=0.0225*(-1.*t247*t270 + t340) - 0.0265*t363 + t337*t367 + t270*t371;
  p_output1[9]=-0.0265*t290 + 0.0225*(-1.*t244*t247 + t362) + t270*t367 + t244*t371;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
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

#include "J_rightToe.hh"

namespace SymFunction
{

void J_rightToe_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
