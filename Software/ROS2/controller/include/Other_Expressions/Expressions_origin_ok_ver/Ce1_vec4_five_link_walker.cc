/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:30 GMT-05:00
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
  double t650;
  double t725;
  double t626;
  double t643;
  double t576;
  double t629;
  double t658;
  double t664;
  double t716;
  double t723;
  double t724;
  double t768;
  double t769;
  double t770;
  double t627;
  double t647;
  double t648;
  double t649;
  double t737;
  double t752;
  double t753;
  double t787;
  double t788;
  double t789;
  double t776;
  double t783;
  double t784;
  double t785;
  double t786;
  double t790;
  double t800;
  double t801;
  double t803;
  double t754;
  double t762;
  double t766;
  double t804;
  double t805;
  double t806;
  double t807;
  double t809;
  double t813;
  double t791;
  double t831;
  double t832;
  double t833;
  double t834;
  double t792;
  double t835;
  double t820;
  double t821;
  double t822;
  double t767;
  double t781;
  double t851;
  double t816;
  double t817;
  double t818;
  double t852;
  double t853;
  double t854;
  double t878;
  double t879;
  double t880;
  double t863;
  double t864;
  double t865;
  double t867;
  double t868;
  double t873;
  double t877;
  double t881;
  double t901;
  double t902;
  double t886;
  double t904;
  double t905;
  double t888;
  t650 = Cos(var1[4]);
  t725 = Sin(var1[4]);
  t626 = Sin(var1[2]);
  t643 = Sin(var1[3]);
  t576 = Cos(var1[3]);
  t629 = Cos(var1[2]);
  t658 = -1.*t650;
  t664 = 1. + t658;
  t716 = 0.5*t664;
  t723 = 0.671885*t650;
  t724 = t716 + t723;
  t768 = t576*t650;
  t769 = -1.*t643*t725;
  t770 = t768 + t769;
  t627 = -1.*t576*t626;
  t647 = -1.*t629*t643;
  t648 = t627 + t647;
  t649 = 0.51185934*t648;
  t737 = t724*t725;
  t752 = -0.171885*t650*t725;
  t753 = t737 + t752;
  t787 = -1.*t650*t643;
  t788 = -1.*t576*t725;
  t789 = t787 + t788;
  t776 = t629*t770;
  t783 = t724*t650;
  t784 = Power(t725,2);
  t785 = 0.171885*t784;
  t786 = t783 + t785;
  t790 = t629*t789;
  t800 = t626*t789;
  t801 = t800 + t776;
  t803 = 0.85216*t753*t801;
  t754 = t650*t643;
  t762 = t576*t725;
  t766 = t754 + t762;
  t804 = -1.*t576*t650;
  t805 = t643*t725;
  t806 = t804 + t805;
  t807 = t626*t806;
  t809 = t790 + t807;
  t813 = 0.85216*t786*t809;
  t791 = -1.*t626*t770;
  t831 = -1.*t629*t576;
  t832 = t626*t643;
  t833 = t831 + t832;
  t834 = 0.51185934*t833;
  t792 = t790 + t791;
  t835 = -1.*t626*t789;
  t820 = Power(t650,2);
  t821 = -0.171885*t820;
  t822 = t783 + t821;
  t767 = -1.*t626*t766;
  t781 = t767 + t776;
  t851 = 0.85216*t753*t792;
  t816 = -1.*t724*t725;
  t817 = 0.171885*t650*t725;
  t818 = t816 + t817;
  t852 = t629*t806;
  t853 = t835 + t852;
  t854 = 0.85216*t786*t853;
  t878 = t576*t724;
  t879 = -0.171885*t643*t725;
  t880 = t878 + t879;
  t863 = -1.*t724*t643;
  t864 = -0.171885*t576*t725;
  t865 = t863 + t864;
  t867 = t724*t643;
  t868 = 0.171885*t576*t725;
  t873 = t867 + t868;
  t877 = t873*t770;
  t881 = t789*t880;
  t901 = -0.171885*t650*t643;
  t902 = t901 + t864;
  t886 = -1.*t789*t873;
  t904 = 0.171885*t576*t650;
  t905 = t904 + t879;
  t888 = -1.*t880*t806;
  p_output1[0]=var2[3]*(-0.5*(t649 + 0.85216*t753*t781 + 0.85216*t786*t792)*var2[2] - 0.5*(t649 + t803 + t813)*var2[3] - 0.5*(t803 + t813 + 0.85216*t801*t818 + 0.85216*(t629*t766 + t626*t770)*t822)*var2[4]);
  p_output1[1]=var2[3]*(-0.5*(0.85216*t753*(-1.*t629*t766 + t791) + t834 + 0.85216*t786*(-1.*t629*t770 + t835))*var2[2] - 0.5*(t834 + t851 + t854)*var2[3] - 0.5*(0.85216*t792*t818 + 0.85216*t781*t822 + t851 + t854)*var2[4]);
  p_output1[2]=var2[3]*(-0.5*(0.85216*t786*(t770*t865 + t877 + t766*t880 + t881) + 0.85216*t753*(-1.*t789*t865 - 1.*t770*t880 + t886 + t888))*var2[3] - 0.5*(0.85216*t818*(t766*t873 + t770*t880) + 0.85216*t822*(-1.*t770*t873 - 1.*t789*t880) + 0.85216*t786*(t877 + t881 + t770*t902 + t766*t905) + 0.85216*t753*(t886 + t888 - 1.*t789*t902 - 1.*t770*t905))*var2[4]);
  p_output1[3]=-0.5*(1.70432*t786*t818 + 1.70432*t753*t822)*var2[3]*var2[4];
  p_output1[4]=-0.0732367608*t818*var2[3]*var2[4];
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

#include "Ce1_vec4_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
