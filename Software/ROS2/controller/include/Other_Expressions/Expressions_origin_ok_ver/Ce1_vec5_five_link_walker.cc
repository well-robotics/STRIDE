/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:31 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
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


#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t814;
  double t793;
  double t794;
  double t815;
  double t782;
  double t799;
  double t819;
  double t823;
  double t824;
  double t825;
  double t837;
  double t838;
  double t839;
  double t844;
  double t848;
  double t826;
  double t827;
  double t828;
  double t857;
  double t861;
  double t862;
  double t885;
  double t887;
  double t889;
  double t890;
  double t891;
  double t900;
  double t903;
  double t906;
  double t893;
  double t896;
  double t897;
  double t898;
  double t899;
  double t907;
  double t908;
  double t909;
  double t910;
  t814 = Cos(var1[3]);
  t793 = Cos(var1[4]);
  t794 = Sin(var1[3]);
  t815 = Sin(var1[4]);
  t782 = Cos(var1[2]);
  t799 = -1.*t793*t794;
  t819 = -1.*t814*t815;
  t823 = t799 + t819;
  t824 = t782*t823;
  t825 = Sin(var1[2]);
  t837 = -1.*t814*t793;
  t838 = t794*t815;
  t839 = t837 + t838;
  t844 = t825*t839;
  t848 = t824 + t844;
  t826 = t814*t793;
  t827 = -1.*t794*t815;
  t828 = t826 + t827;
  t857 = -1.*t825*t823;
  t861 = t782*t839;
  t862 = t857 + t861;
  t885 = -1.*t793;
  t887 = 1. + t885;
  t889 = 0.5*t887;
  t890 = 0.671885*t793;
  t891 = t889 + t890;
  t900 = t814*t891;
  t903 = -0.171885*t794*t815;
  t906 = t900 + t903;
  t893 = -0.171885*t814*t815;
  t896 = t891*t794;
  t897 = 0.171885*t814*t815;
  t898 = t896 + t897;
  t899 = t898*t828;
  t907 = t823*t906;
  t908 = t793*t794;
  t909 = t814*t815;
  t910 = t908 + t909;
  p_output1[0]=var2[4]*(-0.0732367608*(t824 - 1.*t825*t828)*var2[2] - 0.0732367608*t848*var2[3] - 0.0732367608*t848*var2[4]);
  p_output1[1]=var2[4]*(-0.0732367608*(-1.*t782*t828 + t857)*var2[2] - 0.0732367608*t862*var2[3] - 0.0732367608*t862*var2[4]);
  p_output1[2]=var2[4]*(-0.0732367608*(t828*(-1.*t794*t891 + t893) + t899 + t907 + t906*t910)*var2[3] - 0.0732367608*(t828*(-0.171885*t793*t794 + t893) + t899 + t907 + (0.171885*t793*t814 + t903)*t910)*var2[4]);
  p_output1[3]=-0.0732367608*(0.171885*t793*t815 - 1.*t815*t891)*Power(var2[4],2);
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

#include "Ce1_vec5_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
