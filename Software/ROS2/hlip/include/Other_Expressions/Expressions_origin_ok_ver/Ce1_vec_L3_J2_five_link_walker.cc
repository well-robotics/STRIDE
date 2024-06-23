/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:24 GMT-05:00
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
  double t163;
  double t164;
  double t168;
  double t172;
  double t175;
  double t178;
  double t181;
  double t191;
  double t193;
  double t194;
  double t182;
  double t185;
  double t187;
  double t188;
  double t189;
  double t190;
  double t197;
  double t200;
  double t201;
  double t203;
  double t205;
  double t206;
  double t207;
  double t212;
  double t213;
  double t214;
  double t227;
  double t228;
  double t229;
  double t223;
  double t224;
  double t225;
  double t226;
  double t230;
  double t231;
  double t232;
  double t233;
  double t234;
  double t235;
  double t236;
  double t240;
  double t242;
  double t243;
  double t254;
  double t255;
  double t256;
  t163 = Cos(var1[5]);
  t164 = Sin(var1[2]);
  t168 = -1.*t163*t164;
  t172 = Cos(var1[2]);
  t175 = Sin(var1[5]);
  t178 = -1.*t172*t175;
  t181 = t168 + t178;
  t191 = t172*t163;
  t193 = -1.*t164*t175;
  t194 = t191 + t193;
  t182 = Power(t181,2);
  t185 = 0.64429*t182;
  t187 = t163*t164;
  t188 = t172*t175;
  t189 = t187 + t188;
  t190 = 0.64429*t181*t189;
  t197 = Power(t194,2);
  t200 = 0.64429*t197;
  t201 = -1.*t172*t163;
  t203 = t164*t175;
  t205 = t201 + t203;
  t206 = 0.64429*t194*t205;
  t207 = t185 + t190 + t200 + t206;
  t212 = 1.28858*t181*t194;
  t213 = 1.28858*t181*t205;
  t214 = t212 + t213;
  t227 = -0.001112*t163;
  t228 = -0.078865*t175;
  t229 = t227 + t228;
  t223 = -0.078865*t163;
  t224 = 0.001112*t175;
  t225 = t223 + t224;
  t226 = t163*t225;
  t230 = t229*t175;
  t231 = t226 + t230;
  t232 = 0.64429*t205*t231;
  t233 = -1.*t163*t229;
  t234 = t225*t175;
  t235 = t233 + t234;
  t236 = 0.64429*t181*t235;
  t240 = 0.001112*t163;
  t242 = 0.078865*t175;
  t243 = t240 + t242;
  t254 = 0.00071645048*t181;
  t255 = -0.050811930850000006*t205;
  t256 = t254 + t255;
  p_output1[0]=var2[1]*(-0.5*t207*var2[2] - 0.5*t207*var2[5]);
  p_output1[1]=var2[1]*(-0.5*t214*var2[2] - 0.5*t214*var2[5]);
  p_output1[2]=var2[1]*(-0.5*(t232 + t236)*var2[2] - 0.5*(t232 + t236 + 0.64429*t181*(t163*t229 + t163*t243) + 0.64429*t194*(t230 + t175*t243))*var2[5]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[1]*(-0.5*t256*var2[2] - 0.5*t256*var2[5]);
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

#include "Ce1_vec_L3_J2_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L3_J2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
