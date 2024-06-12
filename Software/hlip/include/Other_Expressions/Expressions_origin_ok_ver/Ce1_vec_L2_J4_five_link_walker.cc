/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:13 GMT-05:00
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
  double t118;
  double t107;
  double t111;
  double t124;
  double t115;
  double t126;
  double t133;
  double t139;
  double t142;
  double t144;
  double t147;
  double t153;
  double t157;
  double t168;
  double t170;
  double t172;
  double t175;
  double t178;
  double t179;
  double t188;
  double t189;
  double t190;
  double t193;
  double t194;
  double t197;
  t118 = Cos(var1[2]);
  t107 = Cos(var1[3]);
  t111 = Sin(var1[2]);
  t124 = Sin(var1[3]);
  t115 = -1.*t107*t111;
  t126 = -1.*t118*t124;
  t133 = t115 + t126;
  t139 = -0.05089692188*t133;
  t142 = t118*t107;
  t144 = -1.*t111*t124;
  t147 = t142 + t144;
  t153 = 0.0006594708000000001*t147;
  t157 = t139 + t153;
  t168 = 0.0006594708000000001*t133;
  t170 = -1.*t118*t107;
  t172 = t111*t124;
  t175 = t170 + t172;
  t178 = -0.05089692188*t175;
  t179 = t168 + t178;
  t188 = -0.00102*t107;
  t189 = -0.078722*t124;
  t190 = t188 + t189;
  t193 = 0.00102*t107;
  t194 = 0.078722*t124;
  t197 = t193 + t194;
  p_output1[0]=var2[3]*(-0.5*t157*var2[2] - 0.5*t157*var2[3]);
  p_output1[1]=var2[3]*(-0.5*t179*var2[2] - 0.5*t179*var2[3]);
  p_output1[2]=-0.5*(-0.05089692188*(t107*t190 + t107*t197) + 0.0006594708000000001*(t124*t190 + t124*t197))*Power(var2[3],2);
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

#include "Ce1_vec_L2_J4_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L2_J4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
