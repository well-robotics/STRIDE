/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:39 GMT-05:00
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
  double t105;
  double t120;
  double t133;
  double t136;
  double t112;
  double t93;
  double t146;
  double t149;
  double t151;
  double t152;
  double t138;
  double t142;
  double t143;
  double t144;
  double t145;
  double t153;
  double t155;
  double t111;
  double t124;
  double t125;
  double t169;
  double t170;
  double t174;
  double t210;
  double t211;
  double t212;
  double t203;
  double t204;
  double t205;
  double t197;
  double t198;
  double t159;
  double t160;
  double t168;
  double t184;
  double t190;
  double t191;
  double t199;
  double t207;
  double t208;
  double t215;
  double t216;
  double t218;
  double t220;
  double t221;
  t105 = Cos(var1[3]);
  t120 = Sin(var1[3]);
  t133 = -1.*t105;
  t136 = 1. + t133;
  t112 = Sin(var1[2]);
  t93 = Cos(var1[2]);
  t146 = -0.0695*t136;
  t149 = -0.15232*t105;
  t151 = 0.0011329999999999986*t120;
  t152 = t146 + t149 + t151;
  t138 = -0.0265*t136;
  t142 = -0.025367*t105;
  t143 = 0.08282*t120;
  t144 = t138 + t142 + t143;
  t145 = -1.*t105*t144;
  t153 = -1.*t152*t120;
  t155 = t145 + t153;
  t111 = t93*t105;
  t124 = t112*t120;
  t125 = t111 + t124;
  t169 = t105*t152;
  t170 = -1.*t144*t120;
  t174 = t169 + t170;
  t210 = t105*t112;
  t211 = -1.*t93*t120;
  t212 = t210 + t211;
  t203 = 0.08282*t105;
  t204 = -0.0011329999999999986*t120;
  t205 = t203 + t204;
  t197 = 0.0011329999999999986*t105;
  t198 = t197 + t143;
  t159 = -1.*t105*t112;
  t160 = t93*t120;
  t168 = t159 + t160;
  t184 = -1.*t93*t105;
  t190 = -1.*t112*t120;
  t191 = t184 + t190;
  t199 = t105*t198;
  t207 = -1.*t205*t120;
  t208 = t145 + t199 + t207 + t153;
  t215 = -1.*t105*t205;
  t216 = -1.*t105*t152;
  t218 = t144*t120;
  t220 = -1.*t198*t120;
  t221 = t215 + t216 + t218 + t220;
  p_output1[0]=var2[2]*(-0.5*(0.69051*t125*t155 + 0.69051*t168*t174)*var2[2] - 0.5*(0.69051*t155*t191 + 0.69051*t125*t208 + 0.69051*t174*t212 + 0.69051*t212*t221)*var2[3]);
  p_output1[1]=var2[2]*(-0.5*(0.69051*t155*t168 + 0.69051*t174*t191)*var2[2] - 0.5*(0.69051*t125*t174 + 0.69051*t168*t208 + 0.69051*t155*t212 + 0.69051*t125*t221)*var2[3]);
  p_output1[2]=-0.5*(1.38102*t174*t208 + 1.38102*t155*t221)*var2[2]*var2[3];
  p_output1[3]=-0.5*(0.0571880382*t208 + 0.0007823478299999989*t221)*var2[2]*var2[3];
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

#include "Ce1_vec_L2_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L2_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
