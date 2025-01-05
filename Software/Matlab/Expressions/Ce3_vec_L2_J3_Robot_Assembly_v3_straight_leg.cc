/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:15 GMT-05:00
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
  double t1020;
  double t1042;
  double t1050;
  double t1053;
  double t1040;
  double t1018;
  double t1066;
  double t1076;
  double t1079;
  double t1082;
  double t1054;
  double t1056;
  double t1061;
  double t1063;
  double t1093;
  double t1097;
  double t1098;
  double t1064;
  double t1084;
  double t1088;
  double t1108;
  double t1109;
  double t1110;
  double t1150;
  double t1151;
  double t1153;
  double t1145;
  double t1147;
  double t1124;
  double t1126;
  double t1131;
  double t1021;
  double t1044;
  double t1045;
  double t1148;
  double t1155;
  double t1156;
  double t1176;
  double t1177;
  double t1178;
  double t1158;
  double t1159;
  double t1160;
  double t1161;
  double t1167;
  t1020 = Cos(var1[3]);
  t1042 = Sin(var1[3]);
  t1050 = -1.*t1020;
  t1053 = 1. + t1050;
  t1040 = Sin(var1[2]);
  t1018 = Cos(var1[2]);
  t1066 = -0.0695*t1053;
  t1076 = -0.15232*t1020;
  t1079 = 0.0011329999999999986*t1042;
  t1082 = t1066 + t1076 + t1079;
  t1054 = -0.0265*t1053;
  t1056 = -0.025367*t1020;
  t1061 = 0.08282*t1042;
  t1063 = t1054 + t1056 + t1061;
  t1093 = -1.*t1020*t1040;
  t1097 = t1018*t1042;
  t1098 = t1093 + t1097;
  t1064 = -1.*t1020*t1063;
  t1084 = -1.*t1082*t1042;
  t1088 = t1064 + t1084;
  t1108 = t1020*t1082;
  t1109 = -1.*t1063*t1042;
  t1110 = t1108 + t1109;
  t1150 = 0.08282*t1020;
  t1151 = -0.0011329999999999986*t1042;
  t1153 = t1150 + t1151;
  t1145 = 0.0011329999999999986*t1020;
  t1147 = t1145 + t1061;
  t1124 = -1.*t1018*t1020;
  t1126 = -1.*t1040*t1042;
  t1131 = t1124 + t1126;
  t1021 = t1018*t1020;
  t1044 = t1040*t1042;
  t1045 = t1021 + t1044;
  t1148 = t1020*t1147;
  t1155 = -1.*t1153*t1042;
  t1156 = t1064 + t1148 + t1155 + t1084;
  t1176 = t1020*t1040;
  t1177 = -1.*t1018*t1042;
  t1178 = t1176 + t1177;
  t1158 = -1.*t1020*t1153;
  t1159 = -1.*t1020*t1082;
  t1160 = t1063*t1042;
  t1161 = -1.*t1147*t1042;
  t1167 = t1158 + t1159 + t1160 + t1161;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.69051*t1045*t1088 + 0.69051*t1098*t1110)*var2[0] + 0.5*(0.69051*t1088*t1098 + 0.69051*t1110*t1131)*var2[1])*var2[2];
  p_output1[3]=var2[2]*(0.5*(0.69051*t1088*t1131 + 0.69051*t1045*t1156 + 0.69051*t1110*t1178 + 0.69051*t1167*t1178)*var2[0] + 0.5*(0.69051*t1045*t1110 + 0.69051*t1098*t1156 + 0.69051*t1045*t1167 + 0.69051*t1088*t1178)*var2[1] + 0.5*(1.38102*t1110*t1156 + 1.38102*t1088*t1167)*var2[2] + 0.5*(0.0571880382*t1156 + 0.0007823478299999989*t1167)*var2[3]);
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

#include "Ce3_vec_L2_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L2_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
