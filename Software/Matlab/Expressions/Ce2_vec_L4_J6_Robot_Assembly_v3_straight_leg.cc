/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:32:07 GMT-05:00
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
  double t856;
  double t845;
  double t847;
  double t858;
  double t851;
  double t864;
  double t868;
  double t874;
  double t876;
  double t881;
  double t888;
  double t890;
  double t892;
  double t935;
  double t937;
  double t948;
  double t949;
  double t953;
  double t954;
  double t939;
  double t941;
  double t943;
  double t944;
  double t947;
  double t980;
  double t981;
  double t956;
  double t973;
  double t974;
  double t976;
  double t925;
  double t926;
  double t927;
  double t887;
  double t909;
  double t918;
  double t920;
  double t921;
  double t924;
  double t928;
  double t933;
  double t957;
  double t963;
  double t964;
  double t969;
  double t972;
  double t977;
  double t978;
  double t985;
  double t989;
  double t993;
  double t994;
  double t995;
  double t959;
  double t970;
  double t990;
  double t996;
  double t997;
  double t1011;
  double t1012;
  double t1014;
  double t1015;
  double t1016;
  double t869;
  double t882;
  double t883;
  double t1002;
  double t1003;
  double t1004;
  double t1020;
  double t1021;
  double t1022;
  t856 = Cos(var1[2]);
  t845 = Cos(var1[5]);
  t847 = Sin(var1[2]);
  t858 = Sin(var1[5]);
  t851 = -1.*t845*t847;
  t864 = -1.*t856*t858;
  t868 = t851 + t864;
  t874 = t856*t845;
  t876 = -1.*t847*t858;
  t881 = t874 + t876;
  t888 = t845*t847;
  t890 = t856*t858;
  t892 = t888 + t890;
  t935 = -1.*t845;
  t937 = 1. + t935;
  t948 = -0.0265*t937;
  t949 = -0.025413*t845;
  t953 = -0.08282*t858;
  t954 = t948 + t949 + t953;
  t939 = -0.0695*t937;
  t941 = -0.15232*t845;
  t943 = -0.0010869999999999977*t858;
  t944 = t939 + t941 + t943;
  t947 = t845*t944;
  t980 = -0.08282*t845;
  t981 = t980 + t943;
  t956 = t954*t858;
  t973 = -0.0010869999999999977*t845;
  t974 = 0.08282*t858;
  t976 = t973 + t974;
  t925 = -1.*t856*t845;
  t926 = t847*t858;
  t927 = t925 + t926;
  t887 = 1.38102*t868*t881;
  t909 = Power(t868,2);
  t918 = 0.69051*t909;
  t920 = 0.69051*t868*t892;
  t921 = Power(t881,2);
  t924 = 0.69051*t921;
  t928 = 0.69051*t881*t927;
  t933 = t918 + t920 + t924 + t928;
  t957 = t947 + t956;
  t963 = -1.*t845*t954;
  t964 = t944*t858;
  t969 = t963 + t964;
  t972 = t845*t954;
  t977 = t845*t976;
  t978 = -1.*t944*t858;
  t985 = t981*t858;
  t989 = t972 + t977 + t978 + t985;
  t993 = -1.*t845*t981;
  t994 = t976*t858;
  t995 = t947 + t993 + t956 + t994;
  t959 = 0.69051*t868*t957;
  t970 = 0.69051*t881*t969;
  t990 = 0.69051*t881*t989;
  t996 = 0.69051*t892*t995;
  t997 = t959 + t970 + t990 + t996;
  t1011 = 0.69051*t927*t957;
  t1012 = 0.69051*t868*t969;
  t1014 = 0.69051*t868*t989;
  t1015 = 0.69051*t881*t995;
  t1016 = t1011 + t1012 + t1014 + t1015;
  t869 = -0.0571880382*t868;
  t882 = -0.0007505843699999984*t881;
  t883 = t869 + t882;
  t1002 = -0.0007505843699999984*t868;
  t1003 = -0.0571880382*t927;
  t1004 = t1002 + t1003;
  t1020 = -0.0571880382*t989;
  t1021 = -0.0007505843699999984*t995;
  t1022 = t1020 + t1021;
  p_output1[0]=var2[5]*(-0.5*(t887 + 1.38102*t881*t892)*var2[0] - 0.5*t933*var2[1] - 0.5*t997*var2[2] - 0.5*t883*var2[5]);
  p_output1[1]=var2[5]*(-0.5*t933*var2[0] - 0.5*(t887 + 1.38102*t868*t927)*var2[1] - 0.5*t1016*var2[2] - 0.5*t1004*var2[5]);
  p_output1[2]=var2[5]*(-0.5*t997*var2[0] - 0.5*t1016*var2[1] - 0.5*(1.38102*t957*t989 + 1.38102*t969*t995)*var2[2] - 0.5*t1022*var2[5]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t883*var2[0] - 0.5*t1004*var2[1] - 0.5*t1022*var2[2])*var2[5];
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

#include "Ce2_vec_L4_J6_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L4_J6_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
