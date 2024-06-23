/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:27 GMT-05:00
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
  double t181;
  double t188;
  double t187;
  double t182;
  double t203;
  double t205;
  double t206;
  double t191;
  double t194;
  double t197;
  double t211;
  double t212;
  double t213;
  double t185;
  double t189;
  double t190;
  double t200;
  double t208;
  double t209;
  double t210;
  double t216;
  double t217;
  double t220;
  double t221;
  double t226;
  double t227;
  double t228;
  double t225;
  double t230;
  double t231;
  double t247;
  double t248;
  double t249;
  double t250;
  double t251;
  double t239;
  double t240;
  t181 = Cos(var1[5]);
  t188 = Sin(var1[5]);
  t187 = Cos(var1[2]);
  t182 = Sin(var1[2]);
  t203 = -0.001112*t181;
  t205 = -0.078865*t188;
  t206 = t203 + t205;
  t191 = -0.078865*t181;
  t194 = 0.001112*t188;
  t197 = t191 + t194;
  t211 = t187*t181;
  t212 = -1.*t182*t188;
  t213 = t211 + t212;
  t185 = -1.*t181*t182;
  t189 = -1.*t187*t188;
  t190 = t185 + t189;
  t200 = t181*t197;
  t208 = t206*t188;
  t209 = t200 + t208;
  t210 = 0.64429*t190*t209;
  t216 = -1.*t181*t206;
  t217 = t197*t188;
  t220 = t216 + t217;
  t221 = 0.64429*t213*t220;
  t226 = 0.001112*t181;
  t227 = 0.078865*t188;
  t228 = t226 + t227;
  t225 = t181*t206;
  t230 = t181*t228;
  t231 = t225 + t230;
  t247 = -1.*t187*t181;
  t248 = t182*t188;
  t249 = t247 + t248;
  t250 = 0.64429*t249*t209;
  t251 = 0.64429*t190*t220;
  t239 = t228*t188;
  t240 = t208 + t239;
  p_output1[0]=var2[2]*(-0.5*(t210 + t221)*var2[2] - 0.5*(t210 + t221 + 0.64429*t213*t231 + 0.64429*(t181*t182 + t187*t188)*t240)*var2[5]);
  p_output1[1]=var2[2]*(-0.5*(t250 + t251)*var2[2] - 0.5*(0.64429*t190*t231 + 0.64429*t213*t240 + t250 + t251)*var2[5]);
  p_output1[2]=-0.5*(1.28858*t209*t231 + 1.28858*t220*t240)*var2[2]*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(-0.050811930850000006*t231 + 0.00071645048*t240)*var2[2]*var2[5];
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

#include "Ce1_vec_L3_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L3_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
