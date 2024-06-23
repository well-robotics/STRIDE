/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:40:22 GMT-05:00
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
  double t883;
  double t869;
  double t877;
  double t889;
  double t879;
  double t909;
  double t910;
  double t916;
  double t918;
  double t923;
  double t935;
  double t941;
  double t942;
  double t987;
  double t990;
  double t994;
  double t981;
  double t982;
  double t984;
  double t911;
  double t925;
  double t926;
  double t928;
  double t943;
  double t951;
  double t952;
  double t953;
  double t955;
  double t960;
  double t961;
  double t964;
  double t965;
  double t970;
  double t971;
  double t975;
  double t976;
  double t978;
  double t979;
  double t985;
  double t995;
  double t998;
  double t999;
  double t1000;
  double t1001;
  double t1002;
  double t1004;
  double t1012;
  double t1013;
  double t1014;
  t883 = Cos(var1[2]);
  t869 = Cos(var1[3]);
  t877 = Sin(var1[2]);
  t889 = Sin(var1[3]);
  t879 = -1.*t869*t877;
  t909 = -1.*t883*t889;
  t910 = t879 + t909;
  t916 = -1.*t883*t869;
  t918 = t877*t889;
  t923 = t916 + t918;
  t935 = t883*t869;
  t941 = -1.*t877*t889;
  t942 = t935 + t941;
  t987 = -0.00102*t869;
  t990 = -0.078722*t889;
  t994 = t987 + t990;
  t981 = -0.078722*t869;
  t982 = 0.00102*t889;
  t984 = t981 + t982;
  t911 = 0.0006594708000000001*t910;
  t925 = -0.05089692188*t923;
  t926 = t911 + t925;
  t928 = 0.5*var2[3]*t926;
  t943 = 1.29308*t910*t942;
  t951 = 1.29308*t910*t923;
  t952 = t943 + t951;
  t953 = 0.5*var2[1]*t952;
  t955 = Power(t910,2);
  t960 = 0.64654*t955;
  t961 = t869*t877;
  t964 = t883*t889;
  t965 = t961 + t964;
  t970 = 0.64654*t910*t965;
  t971 = Power(t942,2);
  t975 = 0.64654*t971;
  t976 = 0.64654*t942*t923;
  t978 = t960 + t970 + t975 + t976;
  t979 = 0.5*var2[0]*t978;
  t985 = t869*t984;
  t995 = t994*t889;
  t998 = t985 + t995;
  t999 = 0.64654*t923*t998;
  t1000 = -1.*t869*t994;
  t1001 = t984*t889;
  t1002 = t1000 + t1001;
  t1004 = 0.64654*t910*t1002;
  t1012 = 0.00102*t869;
  t1013 = 0.078722*t889;
  t1014 = t1012 + t1013;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(t928 + t953 + t979 + 0.5*(t1004 + t999)*var2[2]);
  p_output1[3]=var2[1]*(t928 + t953 + t979 + 0.5*(t1004 + 0.64654*t910*(t1014*t869 + t869*t994) + 0.64654*t942*(t1014*t889 + t995) + t999)*var2[2]);
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

#include "Ce3_vec_L2_J2_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L2_J2_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
