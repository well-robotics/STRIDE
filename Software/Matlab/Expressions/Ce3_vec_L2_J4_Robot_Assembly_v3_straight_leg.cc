/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:17 GMT-05:00
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
  double t1054;
  double t1021;
  double t1044;
  double t1056;
  double t1045;
  double t1061;
  double t1066;
  double t1079;
  double t1088;
  double t1089;
  double t1141;
  double t1144;
  double t1145;
  double t1115;
  double t1118;
  double t1119;
  double t1176;
  double t1170;
  double t1172;
  double t1184;
  double t1185;
  double t1186;
  double t1188;
  double t1189;
  double t1190;
  double t1192;
  double t1173;
  double t1175;
  double t1178;
  double t1181;
  double t1182;
  t1054 = Cos(var1[2]);
  t1021 = Cos(var1[3]);
  t1044 = Sin(var1[2]);
  t1056 = Sin(var1[3]);
  t1045 = -1.*t1021*t1044;
  t1061 = t1054*t1056;
  t1066 = t1045 + t1061;
  t1079 = -1.*t1054*t1021;
  t1088 = -1.*t1044*t1056;
  t1089 = t1079 + t1088;
  t1141 = t1021*t1044;
  t1144 = -1.*t1054*t1056;
  t1145 = t1141 + t1144;
  t1115 = t1054*t1021;
  t1118 = t1044*t1056;
  t1119 = t1115 + t1118;
  t1176 = 0.08282*t1056;
  t1170 = -1.*t1021;
  t1172 = 1. + t1170;
  t1184 = 0.08282*t1021;
  t1185 = -0.0011329999999999986*t1056;
  t1186 = t1184 + t1185;
  t1188 = -0.0695*t1172;
  t1189 = -0.15232*t1021;
  t1190 = 0.0011329999999999986*t1056;
  t1192 = t1188 + t1189 + t1190;
  t1173 = -0.0265*t1172;
  t1175 = -0.025367*t1021;
  t1178 = t1173 + t1175 + t1176;
  t1181 = 0.0011329999999999986*t1021;
  t1182 = t1181 + t1176;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.0571880382*t1066 + 0.0007823478299999989*t1119)*var2[0] + 0.5*(0.0007823478299999989*t1066 + 0.0571880382*t1089)*var2[1])*var2[3];
  p_output1[3]=(0.5*(0.0007823478299999989*t1089 + 0.0571880382*t1145)*var2[0] + 0.5*(0.0571880382*t1119 + 0.0007823478299999989*t1145)*var2[1] + 0.5*(0.0007823478299999989*(t1056*t1178 - 1.*t1056*t1182 - 1.*t1021*t1186 - 1.*t1021*t1192) + 0.0571880382*(-1.*t1021*t1178 + t1021*t1182 - 1.*t1056*t1186 - 1.*t1056*t1192))*var2[2])*var2[3];
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

#include "Ce3_vec_L2_J4_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L2_J4_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
