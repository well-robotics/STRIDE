/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:44 GMT-05:00
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
  double t1044;
  double t1037;
  double t1040;
  double t1048;
  double t1041;
  double t1049;
  double t1050;
  double t1053;
  double t1058;
  double t1059;
  double t1066;
  double t1068;
  double t1070;
  double t1073;
  double t1074;
  double t1078;
  double t1081;
  double t1083;
  double t1084;
  double t1085;
  double t1086;
  double t1089;
  double t1090;
  double t1091;
  double t1093;
  double t1094;
  double t1095;
  t1044 = Cos(var1[2]);
  t1037 = Cos(var1[5]);
  t1040 = Sin(var1[2]);
  t1048 = Sin(var1[5]);
  t1041 = -1.*t1037*t1040;
  t1049 = -1.*t1044*t1048;
  t1050 = t1041 + t1049;
  t1053 = -0.050811930850000006*t1050;
  t1058 = t1044*t1037;
  t1059 = -1.*t1040*t1048;
  t1066 = t1058 + t1059;
  t1068 = 0.00071645048*t1066;
  t1070 = t1053 + t1068;
  t1073 = 0.5*var2[0]*t1070;
  t1074 = 0.00071645048*t1050;
  t1078 = -1.*t1044*t1037;
  t1081 = t1040*t1048;
  t1083 = t1078 + t1081;
  t1084 = -0.050811930850000006*t1083;
  t1085 = t1074 + t1084;
  t1086 = 0.5*var2[1]*t1085;
  t1089 = -0.001112*t1037;
  t1090 = -0.078865*t1048;
  t1091 = t1089 + t1090;
  t1093 = 0.001112*t1037;
  t1094 = 0.078865*t1048;
  t1095 = t1093 + t1094;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(t1073 + t1086)*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(t1073 + t1086 + 0.5*(-0.050811930850000006*(t1037*t1091 + t1037*t1095) + 0.00071645048*(t1048*t1091 + t1048*t1095))*var2[2])*var2[5];
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

#include "Ce3_vec_L3_J6_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L3_J6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
