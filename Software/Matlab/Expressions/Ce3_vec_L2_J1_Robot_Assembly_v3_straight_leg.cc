/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:10 GMT-05:00
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
  double t999;
  double t977;
  double t980;
  double t1000;
  double t1015;
  double t1017;
  double t1018;
  double t985;
  double t1002;
  double t1006;
  double t1033;
  double t1034;
  double t1040;
  double t1077;
  double t1078;
  double t1088;
  double t1089;
  double t1093;
  double t1095;
  double t1079;
  double t1082;
  double t1084;
  double t1085;
  double t1057;
  double t1060;
  double t1061;
  double t1042;
  double t1053;
  double t1063;
  double t1064;
  double t1066;
  double t1087;
  double t1096;
  double t1097;
  double t1099;
  double t1100;
  double t1106;
  double t1144;
  double t1145;
  double t1146;
  double t1139;
  double t1140;
  t999 = Cos(var1[2]);
  t977 = Cos(var1[3]);
  t980 = Sin(var1[2]);
  t1000 = Sin(var1[3]);
  t1015 = t999*t977;
  t1017 = t980*t1000;
  t1018 = t1015 + t1017;
  t985 = -1.*t977*t980;
  t1002 = t999*t1000;
  t1006 = t985 + t1002;
  t1033 = t977*t980;
  t1034 = -1.*t999*t1000;
  t1040 = t1033 + t1034;
  t1077 = -1.*t977;
  t1078 = 1. + t1077;
  t1088 = -0.0695*t1078;
  t1089 = -0.15232*t977;
  t1093 = 0.0011329999999999986*t1000;
  t1095 = t1088 + t1089 + t1093;
  t1079 = -0.0265*t1078;
  t1082 = -0.025367*t977;
  t1084 = 0.08282*t1000;
  t1085 = t1079 + t1082 + t1084;
  t1057 = -1.*t999*t977;
  t1060 = -1.*t980*t1000;
  t1061 = t1057 + t1060;
  t1042 = 1.38102*t1040*t1018;
  t1053 = 0.69051*t1040*t1006;
  t1063 = 0.69051*t1061*t1018;
  t1064 = Power(t1018,2);
  t1066 = 0.69051*t1064;
  t1087 = -1.*t977*t1085;
  t1096 = -1.*t1095*t1000;
  t1097 = t1087 + t1096;
  t1099 = t977*t1095;
  t1100 = -1.*t1085*t1000;
  t1106 = t1099 + t1100;
  t1144 = 0.08282*t977;
  t1145 = -0.0011329999999999986*t1000;
  t1146 = t1144 + t1145;
  t1139 = 0.0011329999999999986*t977;
  t1140 = t1139 + t1084;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(0.5*(1.38102*t1006*t1018 + t1042)*var2[0] + 0.5*(0.69051*Power(t1006,2) + t1053 + t1063 + t1066)*var2[1] + 0.5*(0.69051*t1018*t1097 + 0.69051*t1006*t1106)*var2[2] + 0.5*(0.0571880382*t1006 + 0.0007823478299999989*t1018)*var2[3]);
  p_output1[3]=var2[0]*(0.5*(t1042 + 1.38102*t1040*t1061)*var2[0] + 0.5*(0.69051*Power(t1040,2) + t1053 + t1063 + t1066)*var2[1] + 0.5*(0.69051*t1061*t1097 + 0.69051*t1040*t1106 + 0.69051*t1018*(t1087 + t1096 - 1.*t1000*t1146 + t1140*t977) + 0.69051*t1040*(t1000*t1085 - 1.*t1000*t1140 - 1.*t1095*t977 - 1.*t1146*t977))*var2[2] + 0.5*(0.0571880382*t1040 + 0.0007823478299999989*t1061)*var2[3]);
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

#include "Ce3_vec_L2_J1_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L2_J1_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
