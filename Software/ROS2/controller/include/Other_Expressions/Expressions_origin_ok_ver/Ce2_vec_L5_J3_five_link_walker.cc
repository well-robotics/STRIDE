/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:41 GMT-05:00
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
  double t761;
  double t748;
  double t749;
  double t762;
  double t773;
  double t747;
  double t775;
  double t779;
  double t781;
  double t755;
  double t765;
  double t769;
  double t771;
  double t783;
  double t793;
  double t823;
  double t828;
  double t829;
  double t830;
  double t831;
  double t832;
  double t816;
  double t820;
  double t821;
  double t795;
  double t799;
  double t801;
  double t803;
  double t805;
  double t806;
  double t850;
  double t852;
  double t855;
  double t857;
  double t858;
  double t892;
  double t893;
  double t894;
  double t888;
  double t889;
  double t890;
  double t864;
  double t866;
  double t869;
  double t880;
  double t881;
  double t842;
  double t843;
  double t844;
  double t822;
  double t834;
  double t839;
  double t873;
  double t874;
  double t877;
  double t878;
  double t879;
  double t883;
  double t884;
  double t891;
  double t899;
  double t900;
  double t903;
  double t904;
  double t905;
  double t902;
  double t906;
  double t908;
  double t927;
  double t928;
  double t929;
  double t840;
  double t846;
  double t848;
  double t917;
  double t918;
  double t919;
  double t794;
  double t807;
  double t808;
  double t912;
  double t914;
  double t915;
  t761 = Cos(var1[5]);
  t748 = Cos(var1[6]);
  t749 = Sin(var1[5]);
  t762 = Sin(var1[6]);
  t773 = Cos(var1[2]);
  t747 = Sin(var1[2]);
  t775 = t761*t748;
  t779 = -1.*t749*t762;
  t781 = t775 + t779;
  t755 = t748*t749;
  t765 = t761*t762;
  t769 = t755 + t765;
  t771 = -1.*t747*t769;
  t783 = t773*t781;
  t793 = t771 + t783;
  t823 = -1.*t748;
  t828 = 1. + t823;
  t829 = -0.16*t828;
  t830 = -0.167368*t748;
  t831 = 0.022659*t762;
  t832 = t829 + t830 + t831;
  t816 = -0.022659*t748;
  t820 = -0.007367999999999986*t762;
  t821 = t816 + t820;
  t795 = -1.*t748*t749;
  t799 = -1.*t761*t762;
  t801 = t795 + t799;
  t803 = t773*t801;
  t805 = -1.*t747*t781;
  t806 = t803 + t805;
  t850 = t747*t801;
  t852 = t850 + t783;
  t855 = t773*t769;
  t857 = t747*t781;
  t858 = t855 + t857;
  t892 = -1.*t749*t821;
  t893 = t761*t832;
  t894 = t892 + t893;
  t888 = t761*t821;
  t889 = t749*t832;
  t890 = t888 + t889;
  t864 = -1.*t747*t801;
  t866 = -1.*t773*t781;
  t869 = t864 + t866;
  t880 = -1.*t773*t769;
  t881 = t880 + t805;
  t842 = t748*t832;
  t843 = t821*t762;
  t844 = t842 + t843;
  t822 = -1.*t748*t821;
  t834 = t832*t762;
  t839 = t822 + t834;
  t873 = 0.14994*t869*t852;
  t874 = Power(t793,2);
  t877 = 0.14994*t874;
  t878 = Power(t806,2);
  t879 = 0.14994*t878;
  t883 = 0.14994*t881*t858;
  t884 = t873 + t877 + t879 + t883;
  t891 = t890*t769;
  t899 = t894*t781;
  t900 = t891 + t899;
  t903 = -1.*t894*t801;
  t904 = -1.*t890*t781;
  t905 = t903 + t904;
  t902 = 0.14994*t806*t900;
  t906 = 0.14994*t793*t905;
  t908 = t902 + t906;
  t927 = 0.14994*t869*t900;
  t928 = 0.14994*t881*t905;
  t929 = t927 + t928;
  t840 = 0.14994*t839*t793;
  t846 = 0.14994*t844*t806;
  t848 = t840 + t846;
  t917 = 0.14994*t844*t869;
  t918 = 0.14994*t839*t881;
  t919 = t917 + t918;
  t794 = 0.0033974904599999994*t793;
  t807 = -0.0011047579199999977*t806;
  t808 = t794 + t807;
  t912 = -0.0011047579199999977*t869;
  t914 = 0.0033974904599999994*t881;
  t915 = t912 + t914;
  p_output1[0]=var2[2]*(-0.5*(0.29988*t806*t852 + 0.29988*t793*t858)*var2[0] - 0.5*t884*var2[1] - 0.5*t908*var2[2] - 0.5*t848*var2[5] - 0.5*t808*var2[6]);
  p_output1[1]=var2[2]*(-0.5*t884*var2[0] - 0.5*(0.29988*t806*t869 + 0.29988*t793*t881)*var2[1] - 0.5*t929*var2[2] - 0.5*t919*var2[5] - 0.5*t915*var2[6]);
  p_output1[2]=(-0.5*t908*var2[0] - 0.5*t929*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t848*var2[0] - 0.5*t919*var2[1])*var2[2];
  p_output1[6]=(-0.5*t808*var2[0] - 0.5*t915*var2[1])*var2[2];
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

#include "Ce2_vec_L5_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L5_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
