/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:42 GMT-05:00
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
  double t138;
  double t111;
  double t124;
  double t142;
  double t174;
  double t175;
  double t176;
  double t151;
  double t155;
  double t158;
  double t125;
  double t143;
  double t146;
  double t191;
  double t192;
  double t197;
  double t232;
  double t228;
  double t229;
  double t238;
  double t239;
  double t240;
  double t242;
  double t243;
  double t244;
  double t245;
  double t230;
  double t231;
  double t233;
  double t235;
  double t236;
  t138 = Cos(var1[2]);
  t111 = Cos(var1[3]);
  t124 = Sin(var1[2]);
  t142 = Sin(var1[3]);
  t174 = -1.*t111*t124;
  t175 = t138*t142;
  t176 = t174 + t175;
  t151 = -1.*t138*t111;
  t155 = -1.*t124*t142;
  t158 = t151 + t155;
  t125 = t111*t124;
  t143 = -1.*t138*t142;
  t146 = t125 + t143;
  t191 = t138*t111;
  t192 = t124*t142;
  t197 = t191 + t192;
  t232 = 0.08282*t142;
  t228 = -1.*t111;
  t229 = 1. + t228;
  t238 = 0.08282*t111;
  t239 = -0.0011329999999999986*t142;
  t240 = t238 + t239;
  t242 = -0.0695*t229;
  t243 = -0.15232*t111;
  t244 = 0.0011329999999999986*t142;
  t245 = t242 + t243 + t244;
  t230 = -0.0265*t229;
  t231 = -0.025367*t111;
  t233 = t230 + t231 + t232;
  t235 = 0.0011329999999999986*t111;
  t236 = t235 + t232;
  p_output1[0]=var2[3]*(-0.5*(0.0571880382*t176 + 0.0007823478299999989*t197)*var2[2] - 0.5*(0.0571880382*t146 + 0.0007823478299999989*t158)*var2[3]);
  p_output1[1]=var2[3]*(-0.5*(0.0571880382*t158 + 0.0007823478299999989*t176)*var2[2] - 0.5*(0.0007823478299999989*t146 + 0.0571880382*t197)*var2[3]);
  p_output1[2]=-0.5*(0.0007823478299999989*(t142*t233 - 1.*t142*t236 - 1.*t111*t240 - 1.*t111*t245) + 0.0571880382*(-1.*t111*t233 + t111*t236 - 1.*t142*t240 - 1.*t142*t245))*Power(var2[3],2);
  p_output1[3]=0;
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

#include "Ce1_vec_L2_J4_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L2_J4_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
