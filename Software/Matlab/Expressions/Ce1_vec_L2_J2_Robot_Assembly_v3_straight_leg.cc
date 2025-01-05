/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:37 GMT-05:00
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
  double t70;
  double t86;
  double t93;
  double t95;
  double t100;
  double t101;
  double t105;
  double t136;
  double t138;
  double t142;
  double t116;
  double t120;
  double t124;
  double t125;
  double t129;
  double t130;
  double t133;
  double t143;
  double t144;
  double t145;
  double t168;
  double t185;
  double t189;
  double t197;
  double t198;
  double t199;
  double t201;
  double t190;
  double t191;
  double t192;
  double t193;
  double t194;
  double t202;
  double t203;
  double t205;
  double t206;
  double t207;
  double t215;
  double t216;
  double t217;
  double t212;
  double t213;
  t70 = Cos(var1[3]);
  t86 = Sin(var1[2]);
  t93 = t70*t86;
  t95 = Cos(var1[2]);
  t100 = Sin(var1[3]);
  t101 = -1.*t95*t100;
  t105 = t93 + t101;
  t136 = t95*t70;
  t138 = t86*t100;
  t142 = t136 + t138;
  t116 = -1.*t70*t86;
  t120 = t95*t100;
  t124 = t116 + t120;
  t125 = 0.69051*t105*t124;
  t129 = -1.*t95*t70;
  t130 = -1.*t86*t100;
  t133 = t129 + t130;
  t143 = 0.69051*t133*t142;
  t144 = Power(t142,2);
  t145 = 0.69051*t144;
  t168 = 1.38102*t124*t142;
  t185 = -1.*t70;
  t189 = 1. + t185;
  t197 = -0.0695*t189;
  t198 = -0.15232*t70;
  t199 = 0.0011329999999999986*t100;
  t201 = t197 + t198 + t199;
  t190 = -0.0265*t189;
  t191 = -0.025367*t70;
  t192 = 0.08282*t100;
  t193 = t190 + t191 + t192;
  t194 = -1.*t70*t193;
  t202 = -1.*t201*t100;
  t203 = t194 + t202;
  t205 = t70*t201;
  t206 = -1.*t193*t100;
  t207 = t205 + t206;
  t215 = 0.08282*t70;
  t216 = -0.0011329999999999986*t100;
  t217 = t215 + t216;
  t212 = 0.0011329999999999986*t70;
  t213 = t212 + t192;
  p_output1[0]=var2[1]*(-0.5*(0.69051*Power(t124,2) + t125 + t143 + t145)*var2[2] - 0.5*(0.69051*Power(t105,2) + t125 + t143 + t145)*var2[3]);
  p_output1[1]=var2[1]*(-0.5*(1.38102*t124*t133 + t168)*var2[2] - 0.5*(1.38102*t105*t142 + t168)*var2[3]);
  p_output1[2]=var2[1]*(-0.5*(0.69051*t124*t203 + 0.69051*t133*t207)*var2[2] - 0.5*(0.69051*t105*t203 + 0.69051*t142*t207 + 0.69051*t124*(t194 + t202 - 1.*t100*t217 + t213*t70) + 0.69051*t142*(t100*t193 - 1.*t100*t213 - 1.*t201*t70 - 1.*t217*t70))*var2[3]);
  p_output1[3]=var2[1]*(-0.5*(0.0007823478299999989*t124 + 0.0571880382*t133)*var2[2] - 0.5*(0.0007823478299999989*t105 + 0.0571880382*t142)*var2[3]);
  p_output1[4]=0;
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

#include "Ce1_vec_L2_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L2_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
