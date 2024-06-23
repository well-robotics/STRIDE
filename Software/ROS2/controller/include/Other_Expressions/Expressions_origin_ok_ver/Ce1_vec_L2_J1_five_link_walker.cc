/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:05 GMT-05:00
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
  double t71;
  double t60;
  double t65;
  double t72;
  double t78;
  double t80;
  double t83;
  double t66;
  double t73;
  double t75;
  double t84;
  double t86;
  double t90;
  double t93;
  double t95;
  double t96;
  double t107;
  double t111;
  double t112;
  double t115;
  double t118;
  double t123;
  double t124;
  double t125;
  double t126;
  double t129;
  double t147;
  double t148;
  double t149;
  double t142;
  double t143;
  double t144;
  double t145;
  double t152;
  double t153;
  double t154;
  double t155;
  double t157;
  double t160;
  double t161;
  double t166;
  double t168;
  double t170;
  double t183;
  double t185;
  double t186;
  t71 = Cos(var1[2]);
  t60 = Cos(var1[3]);
  t65 = Sin(var1[2]);
  t72 = Sin(var1[3]);
  t78 = t71*t60;
  t80 = -1.*t65*t72;
  t83 = t78 + t80;
  t66 = -1.*t60*t65;
  t73 = -1.*t71*t72;
  t75 = t66 + t73;
  t84 = 1.29308*t75*t83;
  t86 = t60*t65;
  t90 = t71*t72;
  t93 = t86 + t90;
  t95 = 1.29308*t93*t83;
  t96 = t84 + t95;
  t107 = Power(t75,2);
  t111 = 0.64654*t107;
  t112 = 0.64654*t75*t93;
  t115 = Power(t83,2);
  t118 = 0.64654*t115;
  t123 = -1.*t71*t60;
  t124 = t65*t72;
  t125 = t123 + t124;
  t126 = 0.64654*t83*t125;
  t129 = t111 + t112 + t118 + t126;
  t147 = -0.00102*t60;
  t148 = -0.078722*t72;
  t149 = t147 + t148;
  t142 = -0.078722*t60;
  t143 = 0.00102*t72;
  t144 = t142 + t143;
  t145 = t60*t144;
  t152 = t149*t72;
  t153 = t145 + t152;
  t154 = 0.64654*t75*t153;
  t155 = -1.*t60*t149;
  t157 = t144*t72;
  t160 = t155 + t157;
  t161 = 0.64654*t83*t160;
  t166 = 0.00102*t60;
  t168 = 0.078722*t72;
  t170 = t166 + t168;
  t183 = -0.05089692188*t75;
  t185 = 0.0006594708000000001*t83;
  t186 = t183 + t185;
  p_output1[0]=var2[0]*(-0.5*t96*var2[2] - 0.5*t96*var2[3]);
  p_output1[1]=var2[0]*(-0.5*t129*var2[2] - 0.5*t129*var2[3]);
  p_output1[2]=var2[0]*(-0.5*(t154 + t161)*var2[2] - 0.5*(t154 + t161 + 0.64654*(t149*t60 + t170*t60)*t83 + 0.64654*(t152 + t170*t72)*t93)*var2[3]);
  p_output1[3]=var2[0]*(-0.5*t186*var2[2] - 0.5*t186*var2[3]);
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

#include "Ce1_vec_L2_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L2_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
