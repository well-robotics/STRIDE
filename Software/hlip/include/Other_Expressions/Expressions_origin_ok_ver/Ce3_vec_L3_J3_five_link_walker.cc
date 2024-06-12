/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:38 GMT-05:00
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
  double t1026;
  double t1033;
  double t1032;
  double t1025;
  double t1050;
  double t1053;
  double t1054;
  double t1041;
  double t1044;
  double t1048;
  double t1062;
  double t1065;
  double t1066;
  double t1049;
  double t1058;
  double t1059;
  double t1067;
  double t1068;
  double t1070;
  double t1088;
  double t1089;
  double t1090;
  double t1087;
  double t1091;
  double t1092;
  double t1078;
  double t1080;
  double t1081;
  double t1077;
  double t1082;
  double t1094;
  double t1095;
  double t1029;
  double t1037;
  double t1040;
  double t1061;
  double t1071;
  t1026 = Cos(var1[5]);
  t1033 = Sin(var1[5]);
  t1032 = Sin(var1[2]);
  t1025 = Cos(var1[2]);
  t1050 = -0.001112*t1026;
  t1053 = -0.078865*t1033;
  t1054 = t1050 + t1053;
  t1041 = -0.078865*t1026;
  t1044 = 0.001112*t1033;
  t1048 = t1041 + t1044;
  t1062 = -1.*t1026*t1032;
  t1065 = -1.*t1025*t1033;
  t1066 = t1062 + t1065;
  t1049 = t1026*t1048;
  t1058 = t1054*t1033;
  t1059 = t1049 + t1058;
  t1067 = -1.*t1026*t1054;
  t1068 = t1048*t1033;
  t1070 = t1067 + t1068;
  t1088 = 0.001112*t1026;
  t1089 = 0.078865*t1033;
  t1090 = t1088 + t1089;
  t1087 = t1026*t1054;
  t1091 = t1026*t1090;
  t1092 = t1087 + t1091;
  t1078 = t1025*t1026;
  t1080 = -1.*t1032*t1033;
  t1081 = t1078 + t1080;
  t1077 = 0.64429*t1066*t1059;
  t1082 = 0.64429*t1081*t1070;
  t1094 = t1090*t1033;
  t1095 = t1058 + t1094;
  t1029 = -1.*t1025*t1026;
  t1037 = t1032*t1033;
  t1040 = t1029 + t1037;
  t1061 = 0.64429*t1040*t1059;
  t1071 = 0.64429*t1066*t1070;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(t1077 + t1082)*var2[0] + 0.5*(t1061 + t1071)*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[2]*(0.5*(t1077 + t1082 + 0.64429*t1081*t1092 + 0.64429*(t1026*t1032 + t1025*t1033)*t1095)*var2[0] + 0.5*(t1061 + t1071 + 0.64429*t1066*t1092 + 0.64429*t1081*t1095)*var2[1] + 0.5*(1.28858*t1059*t1092 + 1.28858*t1070*t1095)*var2[2] + 0.5*(-0.050811930850000006*t1092 + 0.00071645048*t1095)*var2[5]);
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

#include "Ce3_vec_L3_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L3_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
