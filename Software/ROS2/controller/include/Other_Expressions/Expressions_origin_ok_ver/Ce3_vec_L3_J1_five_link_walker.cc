/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:34 GMT-05:00
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
  double t982;
  double t961;
  double t964;
  double t985;
  double t970;
  double t998;
  double t1006;
  double t1011;
  double t1014;
  double t1015;
  double t1025;
  double t1026;
  double t1027;
  double t1049;
  double t1050;
  double t1052;
  double t1044;
  double t1046;
  double t1047;
  double t1010;
  double t1020;
  double t1021;
  double t1023;
  double t1024;
  double t1029;
  double t1030;
  double t1031;
  double t1032;
  double t1033;
  double t1035;
  double t1036;
  double t1037;
  double t1038;
  double t1039;
  double t1040;
  double t1041;
  double t1042;
  double t1043;
  double t1048;
  double t1053;
  double t1054;
  double t1055;
  double t1056;
  double t1058;
  double t1059;
  double t1060;
  double t1066;
  double t1067;
  double t1068;
  t982 = Cos(var1[2]);
  t961 = Cos(var1[5]);
  t964 = Sin(var1[2]);
  t985 = Sin(var1[5]);
  t970 = -1.*t961*t964;
  t998 = -1.*t982*t985;
  t1006 = t970 + t998;
  t1011 = t982*t961;
  t1014 = -1.*t964*t985;
  t1015 = t1011 + t1014;
  t1025 = t961*t964;
  t1026 = t982*t985;
  t1027 = t1025 + t1026;
  t1049 = -0.001112*t961;
  t1050 = -0.078865*t985;
  t1052 = t1049 + t1050;
  t1044 = -0.078865*t961;
  t1046 = 0.001112*t985;
  t1047 = t1044 + t1046;
  t1010 = -0.050811930850000006*t1006;
  t1020 = 0.00071645048*t1015;
  t1021 = t1010 + t1020;
  t1023 = 0.5*var2[5]*t1021;
  t1024 = 1.28858*t1006*t1015;
  t1029 = 1.28858*t1027*t1015;
  t1030 = t1024 + t1029;
  t1031 = 0.5*var2[0]*t1030;
  t1032 = Power(t1006,2);
  t1033 = 0.64429*t1032;
  t1035 = 0.64429*t1006*t1027;
  t1036 = Power(t1015,2);
  t1037 = 0.64429*t1036;
  t1038 = -1.*t982*t961;
  t1039 = t964*t985;
  t1040 = t1038 + t1039;
  t1041 = 0.64429*t1015*t1040;
  t1042 = t1033 + t1035 + t1037 + t1041;
  t1043 = 0.5*var2[1]*t1042;
  t1048 = t961*t1047;
  t1053 = t1052*t985;
  t1054 = t1048 + t1053;
  t1055 = 0.64429*t1006*t1054;
  t1056 = -1.*t961*t1052;
  t1058 = t1047*t985;
  t1059 = t1056 + t1058;
  t1060 = 0.64429*t1015*t1059;
  t1066 = 0.001112*t961;
  t1067 = 0.078865*t985;
  t1068 = t1066 + t1067;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(t1023 + t1031 + t1043 + 0.5*(t1055 + t1060)*var2[2]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[0]*(t1023 + t1031 + t1043 + 0.5*(t1055 + t1060 + 0.64429*t1015*(t1052*t961 + t1068*t961) + 0.64429*t1027*(t1053 + t1068*t985))*var2[2]);
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

#include "Ce3_vec_L3_J1_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L3_J1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
