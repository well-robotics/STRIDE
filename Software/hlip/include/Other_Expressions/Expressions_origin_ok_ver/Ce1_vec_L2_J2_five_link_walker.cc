/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:37:07 GMT-05:00
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
  double t75;
  double t78;
  double t83;
  double t84;
  double t86;
  double t90;
  double t93;
  double t111;
  double t112;
  double t115;
  double t95;
  double t99;
  double t101;
  double t105;
  double t106;
  double t107;
  double t118;
  double t124;
  double t125;
  double t126;
  double t133;
  double t135;
  double t138;
  double t145;
  double t147;
  double t148;
  double t168;
  double t170;
  double t171;
  double t162;
  double t163;
  double t164;
  double t166;
  double t172;
  double t175;
  double t176;
  double t177;
  double t178;
  double t179;
  double t180;
  double t185;
  double t187;
  double t188;
  double t200;
  double t201;
  double t202;
  t75 = Cos(var1[3]);
  t78 = Sin(var1[2]);
  t83 = -1.*t75*t78;
  t84 = Cos(var1[2]);
  t86 = Sin(var1[3]);
  t90 = -1.*t84*t86;
  t93 = t83 + t90;
  t111 = t84*t75;
  t112 = -1.*t78*t86;
  t115 = t111 + t112;
  t95 = Power(t93,2);
  t99 = 0.64654*t95;
  t101 = t75*t78;
  t105 = t84*t86;
  t106 = t101 + t105;
  t107 = 0.64654*t93*t106;
  t118 = Power(t115,2);
  t124 = 0.64654*t118;
  t125 = -1.*t84*t75;
  t126 = t78*t86;
  t133 = t125 + t126;
  t135 = 0.64654*t115*t133;
  t138 = t99 + t107 + t124 + t135;
  t145 = 1.29308*t93*t115;
  t147 = 1.29308*t93*t133;
  t148 = t145 + t147;
  t168 = -0.00102*t75;
  t170 = -0.078722*t86;
  t171 = t168 + t170;
  t162 = -0.078722*t75;
  t163 = 0.00102*t86;
  t164 = t162 + t163;
  t166 = t75*t164;
  t172 = t171*t86;
  t175 = t166 + t172;
  t176 = 0.64654*t133*t175;
  t177 = -1.*t75*t171;
  t178 = t164*t86;
  t179 = t177 + t178;
  t180 = 0.64654*t93*t179;
  t185 = 0.00102*t75;
  t187 = 0.078722*t86;
  t188 = t185 + t187;
  t200 = 0.0006594708000000001*t93;
  t201 = -0.05089692188*t133;
  t202 = t200 + t201;
  p_output1[0]=var2[1]*(-0.5*t138*var2[2] - 0.5*t138*var2[3]);
  p_output1[1]=var2[1]*(-0.5*t148*var2[2] - 0.5*t148*var2[3]);
  p_output1[2]=var2[1]*(-0.5*(t176 + t180)*var2[2] - 0.5*(t176 + t180 + 0.64654*t115*(t172 + t188*t86) + 0.64654*(t171*t75 + t188*t75)*t93)*var2[3]);
  p_output1[3]=var2[1]*(-0.5*t202*var2[2] - 0.5*t202*var2[3]);
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

#include "Ce1_vec_L2_J2_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L2_J2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
