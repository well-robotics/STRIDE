/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:34 GMT-05:00
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
  double t197;
  double t190;
  double t191;
  double t200;
  double t194;
  double t203;
  double t205;
  double t208;
  double t209;
  double t211;
  double t213;
  double t217;
  double t220;
  double t227;
  double t228;
  double t230;
  double t231;
  double t234;
  double t235;
  double t243;
  double t244;
  double t245;
  double t248;
  double t249;
  double t252;
  t197 = Cos(var1[2]);
  t190 = Cos(var1[5]);
  t191 = Sin(var1[2]);
  t200 = Sin(var1[5]);
  t194 = -1.*t190*t191;
  t203 = -1.*t197*t200;
  t205 = t194 + t203;
  t208 = -0.050811930850000006*t205;
  t209 = t197*t190;
  t211 = -1.*t191*t200;
  t213 = t209 + t211;
  t217 = 0.00071645048*t213;
  t220 = t208 + t217;
  t227 = 0.00071645048*t205;
  t228 = -1.*t197*t190;
  t230 = t191*t200;
  t231 = t228 + t230;
  t234 = -0.050811930850000006*t231;
  t235 = t227 + t234;
  t243 = -0.001112*t190;
  t244 = -0.078865*t200;
  t245 = t243 + t244;
  t248 = 0.001112*t190;
  t249 = 0.078865*t200;
  t252 = t248 + t249;
  p_output1[0]=var2[5]*(-0.5*t220*var2[2] - 0.5*t220*var2[5]);
  p_output1[1]=var2[5]*(-0.5*t235*var2[2] - 0.5*t235*var2[5]);
  p_output1[2]=-0.5*(-0.050811930850000006*(t190*t245 + t190*t252) + 0.00071645048*(t200*t245 + t200*t252))*Power(var2[5],2);
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

#include "Ce1_vec_L3_J6_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L3_J6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
