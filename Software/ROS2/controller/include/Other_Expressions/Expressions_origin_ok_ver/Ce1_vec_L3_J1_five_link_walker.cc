/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:22 GMT-05:00
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
  double t147;
  double t133;
  double t139;
  double t153;
  double t164;
  double t166;
  double t168;
  double t142;
  double t162;
  double t163;
  double t172;
  double t175;
  double t178;
  double t181;
  double t182;
  double t183;
  double t190;
  double t191;
  double t193;
  double t194;
  double t197;
  double t199;
  double t200;
  double t201;
  double t203;
  double t204;
  double t213;
  double t214;
  double t215;
  double t209;
  double t210;
  double t211;
  double t212;
  double t216;
  double t217;
  double t218;
  double t219;
  double t220;
  double t221;
  double t222;
  double t226;
  double t227;
  double t228;
  double t239;
  double t240;
  double t241;
  t147 = Cos(var1[2]);
  t133 = Cos(var1[5]);
  t139 = Sin(var1[2]);
  t153 = Sin(var1[5]);
  t164 = t147*t133;
  t166 = -1.*t139*t153;
  t168 = t164 + t166;
  t142 = -1.*t133*t139;
  t162 = -1.*t147*t153;
  t163 = t142 + t162;
  t172 = 1.28858*t163*t168;
  t175 = t133*t139;
  t178 = t147*t153;
  t181 = t175 + t178;
  t182 = 1.28858*t181*t168;
  t183 = t172 + t182;
  t190 = Power(t163,2);
  t191 = 0.64429*t190;
  t193 = 0.64429*t163*t181;
  t194 = Power(t168,2);
  t197 = 0.64429*t194;
  t199 = -1.*t147*t133;
  t200 = t139*t153;
  t201 = t199 + t200;
  t203 = 0.64429*t168*t201;
  t204 = t191 + t193 + t197 + t203;
  t213 = -0.001112*t133;
  t214 = -0.078865*t153;
  t215 = t213 + t214;
  t209 = -0.078865*t133;
  t210 = 0.001112*t153;
  t211 = t209 + t210;
  t212 = t133*t211;
  t216 = t215*t153;
  t217 = t212 + t216;
  t218 = 0.64429*t163*t217;
  t219 = -1.*t133*t215;
  t220 = t211*t153;
  t221 = t219 + t220;
  t222 = 0.64429*t168*t221;
  t226 = 0.001112*t133;
  t227 = 0.078865*t153;
  t228 = t226 + t227;
  t239 = -0.050811930850000006*t163;
  t240 = 0.00071645048*t168;
  t241 = t239 + t240;
  p_output1[0]=var2[0]*(-0.5*t183*var2[2] - 0.5*t183*var2[5]);
  p_output1[1]=var2[0]*(-0.5*t204*var2[2] - 0.5*t204*var2[5]);
  p_output1[2]=var2[0]*(-0.5*(t218 + t222)*var2[2] - 0.5*(t218 + t222 + 0.64429*t168*(t133*t215 + t133*t228) + 0.64429*t181*(t216 + t153*t228))*var2[5]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(-0.5*t241*var2[2] - 0.5*t241*var2[5]);
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

#include "Ce1_vec_L3_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L3_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
