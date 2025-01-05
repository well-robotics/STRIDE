/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:29:34 GMT-05:00
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
  double t74;
  double t64;
  double t65;
  double t83;
  double t70;
  double t84;
  double t86;
  double t100;
  double t101;
  double t105;
  double t111;
  double t117;
  double t120;
  double t124;
  double t90;
  double t91;
  double t93;
  double t143;
  double t144;
  double t145;
  double t146;
  double t164;
  double t165;
  double t174;
  double t175;
  double t176;
  double t179;
  double t168;
  double t169;
  double t170;
  double t172;
  double t173;
  double t180;
  double t182;
  double t185;
  double t186;
  double t189;
  double t198;
  double t199;
  double t200;
  double t194;
  double t196;
  t74 = Cos(var1[2]);
  t64 = Cos(var1[3]);
  t65 = Sin(var1[2]);
  t83 = Sin(var1[3]);
  t70 = t64*t65;
  t84 = -1.*t74*t83;
  t86 = t70 + t84;
  t100 = t74*t64;
  t101 = t65*t83;
  t105 = t100 + t101;
  t111 = 1.38102*t86*t105;
  t117 = -1.*t64*t65;
  t120 = t74*t83;
  t124 = t117 + t120;
  t90 = -1.*t74*t64;
  t91 = -1.*t65*t83;
  t93 = t90 + t91;
  t143 = 0.69051*t86*t124;
  t144 = 0.69051*t93*t105;
  t145 = Power(t105,2);
  t146 = 0.69051*t145;
  t164 = -1.*t64;
  t165 = 1. + t164;
  t174 = -0.0695*t165;
  t175 = -0.15232*t64;
  t176 = 0.0011329999999999986*t83;
  t179 = t174 + t175 + t176;
  t168 = -0.0265*t165;
  t169 = -0.025367*t64;
  t170 = 0.08282*t83;
  t172 = t168 + t169 + t170;
  t173 = -1.*t64*t172;
  t180 = -1.*t179*t83;
  t182 = t173 + t180;
  t185 = t64*t179;
  t186 = -1.*t172*t83;
  t189 = t185 + t186;
  t198 = 0.08282*t64;
  t199 = -0.0011329999999999986*t83;
  t200 = t198 + t199;
  t194 = 0.0011329999999999986*t64;
  t196 = t194 + t170;
  p_output1[0]=var2[0]*(-0.5*(t111 + 1.38102*t105*t124)*var2[2] - 0.5*(t111 + 1.38102*t86*t93)*var2[3]);
  p_output1[1]=var2[0]*(-0.5*(0.69051*Power(t124,2) + t143 + t144 + t146)*var2[2] - 0.5*(t143 + t144 + t146 + 0.69051*Power(t86,2))*var2[3]);
  p_output1[2]=var2[0]*(-0.5*(0.69051*t105*t182 + 0.69051*t124*t189)*var2[2] - 0.5*(0.69051*t105*(t173 + t180 + t196*t64 - 1.*t200*t83) + 0.69051*t189*t86 + 0.69051*(-1.*t179*t64 - 1.*t200*t64 + t172*t83 - 1.*t196*t83)*t86 + 0.69051*t182*t93)*var2[3]);
  p_output1[3]=var2[0]*(-0.5*(0.0007823478299999989*t105 + 0.0571880382*t124)*var2[2] - 0.5*(0.0571880382*t86 + 0.0007823478299999989*t93)*var2[3]);
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

#include "Ce1_vec_L2_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L2_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
