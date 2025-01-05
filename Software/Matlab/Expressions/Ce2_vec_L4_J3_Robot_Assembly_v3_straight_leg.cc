/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:32:00 GMT-05:00
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
  double t840;
  double t731;
  double t825;
  double t841;
  double t830;
  double t842;
  double t845;
  double t851;
  double t852;
  double t856;
  double t869;
  double t874;
  double t876;
  double t912;
  double t913;
  double t926;
  double t927;
  double t928;
  double t929;
  double t918;
  double t920;
  double t921;
  double t924;
  double t896;
  double t897;
  double t899;
  double t868;
  double t886;
  double t887;
  double t888;
  double t890;
  double t892;
  double t902;
  double t906;
  double t925;
  double t933;
  double t934;
  double t936;
  double t937;
  double t939;
  double t963;
  double t964;
  double t966;
  double t935;
  double t941;
  double t942;
  double t847;
  double t858;
  double t859;
  double t948;
  double t949;
  double t952;
  t840 = Cos(var1[2]);
  t731 = Cos(var1[5]);
  t825 = Sin(var1[2]);
  t841 = Sin(var1[5]);
  t830 = -1.*t731*t825;
  t842 = -1.*t840*t841;
  t845 = t830 + t842;
  t851 = t840*t731;
  t852 = -1.*t825*t841;
  t856 = t851 + t852;
  t869 = t731*t825;
  t874 = t840*t841;
  t876 = t869 + t874;
  t912 = -1.*t731;
  t913 = 1. + t912;
  t926 = -0.0265*t913;
  t927 = -0.025413*t731;
  t928 = -0.08282*t841;
  t929 = t926 + t927 + t928;
  t918 = -0.0695*t913;
  t920 = -0.15232*t731;
  t921 = -0.0010869999999999977*t841;
  t924 = t918 + t920 + t921;
  t896 = -1.*t840*t731;
  t897 = t825*t841;
  t899 = t896 + t897;
  t868 = 1.38102*t845*t856;
  t886 = Power(t845,2);
  t887 = 0.69051*t886;
  t888 = 0.69051*t845*t876;
  t890 = Power(t856,2);
  t892 = 0.69051*t890;
  t902 = 0.69051*t856*t899;
  t906 = t887 + t888 + t892 + t902;
  t925 = t731*t924;
  t933 = t929*t841;
  t934 = t925 + t933;
  t936 = -1.*t731*t929;
  t937 = t924*t841;
  t939 = t936 + t937;
  t963 = 0.69051*t899*t934;
  t964 = 0.69051*t845*t939;
  t966 = t963 + t964;
  t935 = 0.69051*t845*t934;
  t941 = 0.69051*t856*t939;
  t942 = t935 + t941;
  t847 = -0.0571880382*t845;
  t858 = -0.0007505843699999984*t856;
  t859 = t847 + t858;
  t948 = -0.0007505843699999984*t845;
  t949 = -0.0571880382*t899;
  t952 = t948 + t949;
  p_output1[0]=var2[2]*(-0.5*(t868 + 1.38102*t856*t876)*var2[0] - 0.5*t906*var2[1] - 0.5*t942*var2[2] - 0.5*t859*var2[5]);
  p_output1[1]=var2[2]*(-0.5*t906*var2[0] - 0.5*(t868 + 1.38102*t845*t899)*var2[1] - 0.5*t966*var2[2] - 0.5*t952*var2[5]);
  p_output1[2]=(-0.5*t942*var2[0] - 0.5*t966*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t859*var2[0] - 0.5*t952*var2[1])*var2[2];
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

#include "Ce2_vec_L4_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce2_vec_L4_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
