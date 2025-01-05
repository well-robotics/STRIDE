/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:12 GMT-05:00
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
  double t1015;
  double t1002;
  double t1006;
  double t1017;
  double t1009;
  double t1018;
  double t1020;
  double t1028;
  double t1033;
  double t1040;
  double t1053;
  double t1054;
  double t1056;
  double t1099;
  double t1106;
  double t1115;
  double t1118;
  double t1119;
  double t1121;
  double t1108;
  double t1109;
  double t1110;
  double t1111;
  double t1066;
  double t1070;
  double t1076;
  double t1061;
  double t1079;
  double t1088;
  double t1089;
  double t1093;
  double t1114;
  double t1122;
  double t1123;
  double t1126;
  double t1130;
  double t1131;
  double t1169;
  double t1170;
  double t1171;
  double t1161;
  double t1166;
  t1015 = Cos(var1[2]);
  t1002 = Cos(var1[3]);
  t1006 = Sin(var1[2]);
  t1017 = Sin(var1[3]);
  t1009 = -1.*t1002*t1006;
  t1018 = t1015*t1017;
  t1020 = t1009 + t1018;
  t1028 = -1.*t1015*t1002;
  t1033 = -1.*t1006*t1017;
  t1040 = t1028 + t1033;
  t1053 = t1015*t1002;
  t1054 = t1006*t1017;
  t1056 = t1053 + t1054;
  t1099 = -1.*t1002;
  t1106 = 1. + t1099;
  t1115 = -0.0695*t1106;
  t1118 = -0.15232*t1002;
  t1119 = 0.0011329999999999986*t1017;
  t1121 = t1115 + t1118 + t1119;
  t1108 = -0.0265*t1106;
  t1109 = -0.025367*t1002;
  t1110 = 0.08282*t1017;
  t1111 = t1108 + t1109 + t1110;
  t1066 = t1002*t1006;
  t1070 = -1.*t1015*t1017;
  t1076 = t1066 + t1070;
  t1061 = 1.38102*t1020*t1056;
  t1079 = 0.69051*t1076*t1020;
  t1088 = 0.69051*t1040*t1056;
  t1089 = Power(t1056,2);
  t1093 = 0.69051*t1089;
  t1114 = -1.*t1002*t1111;
  t1122 = -1.*t1121*t1017;
  t1123 = t1114 + t1122;
  t1126 = t1002*t1121;
  t1130 = -1.*t1111*t1017;
  t1131 = t1126 + t1130;
  t1169 = 0.08282*t1002;
  t1170 = -0.0011329999999999986*t1017;
  t1171 = t1169 + t1170;
  t1161 = 0.0011329999999999986*t1002;
  t1166 = t1161 + t1110;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(0.5*(0.69051*Power(t1020,2) + t1079 + t1088 + t1093)*var2[0] + 0.5*(1.38102*t1020*t1040 + t1061)*var2[1] + 0.5*(0.69051*t1020*t1123 + 0.69051*t1040*t1131)*var2[2] + 0.5*(0.0007823478299999989*t1020 + 0.0571880382*t1040)*var2[3]);
  p_output1[3]=var2[1]*(0.5*(0.69051*Power(t1076,2) + t1079 + t1088 + t1093)*var2[0] + 0.5*(t1061 + 1.38102*t1056*t1076)*var2[1] + 0.5*(0.69051*t1076*t1123 + 0.69051*t1056*t1131 + 0.69051*t1056*(t1017*t1111 - 1.*t1002*t1121 - 1.*t1017*t1166 - 1.*t1002*t1171) + 0.69051*t1020*(t1114 + t1122 + t1002*t1166 - 1.*t1017*t1171))*var2[2] + 0.5*(0.0571880382*t1056 + 0.0007823478299999989*t1076)*var2[3]);
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

#include "Ce3_vec_L2_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L2_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
