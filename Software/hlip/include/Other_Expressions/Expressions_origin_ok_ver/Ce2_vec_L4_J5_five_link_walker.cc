/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:30 GMT-05:00
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
  double t738;
  double t726;
  double t728;
  double t739;
  double t748;
  double t736;
  double t741;
  double t744;
  double t722;
  double t747;
  double t749;
  double t753;
  double t754;
  double t755;
  double t761;
  double t781;
  double t783;
  double t786;
  double t790;
  double t791;
  double t793;
  double t794;
  double t795;
  double t797;
  double t816;
  double t817;
  double t806;
  double t807;
  double t808;
  double t823;
  double t828;
  double t765;
  double t767;
  double t769;
  double t770;
  double t771;
  double t773;
  double t831;
  double t832;
  double t833;
  double t834;
  double t839;
  double t840;
  double t856;
  double t857;
  double t867;
  double t869;
  double t870;
  double t873;
  double t874;
  double t875;
  double t885;
  double t886;
  double t887;
  double t880;
  double t881;
  double t882;
  double t859;
  double t860;
  double t861;
  double t827;
  double t829;
  double t830;
  double t853;
  double t854;
  double t787;
  double t799;
  double t801;
  double t805;
  double t810;
  double t813;
  double t820;
  double t821;
  double t843;
  double t855;
  double t858;
  double t863;
  double t864;
  double t865;
  double t871;
  double t876;
  double t877;
  double t879;
  double t883;
  double t884;
  double t888;
  double t889;
  double t891;
  double t892;
  double t893;
  double t895;
  double t896;
  double t897;
  double t898;
  double t899;
  double t920;
  double t921;
  double t922;
  double t923;
  double t924;
  double t878;
  double t890;
  double t894;
  double t900;
  double t901;
  double t909;
  double t910;
  double t911;
  double t912;
  double t913;
  double t803;
  double t822;
  double t842;
  double t844;
  double t845;
  double t932;
  double t933;
  double t934;
  double t935;
  double t936;
  double t946;
  double t947;
  double t948;
  double t905;
  double t906;
  double t907;
  double t762;
  double t775;
  double t777;
  double t928;
  double t929;
  double t930;
  t738 = Cos(var1[3]);
  t726 = Cos(var1[4]);
  t728 = Sin(var1[3]);
  t739 = Sin(var1[4]);
  t748 = Cos(var1[2]);
  t736 = -1.*t726*t728;
  t741 = -1.*t738*t739;
  t744 = t736 + t741;
  t722 = Sin(var1[2]);
  t747 = t722*t744;
  t749 = t738*t726;
  t753 = -1.*t728*t739;
  t754 = t749 + t753;
  t755 = t748*t754;
  t761 = t747 + t755;
  t781 = -0.022663*t726;
  t783 = -0.007370999999999989*t739;
  t786 = t781 + t783;
  t790 = -1.*t726;
  t791 = 1. + t790;
  t793 = -0.16*t791;
  t794 = -0.167371*t726;
  t795 = 0.022663*t739;
  t797 = t793 + t794 + t795;
  t816 = -0.007370999999999989*t726;
  t817 = t816 + t795;
  t806 = 0.022663*t726;
  t807 = 0.007370999999999989*t739;
  t808 = t806 + t807;
  t823 = t726*t797;
  t828 = t786*t739;
  t765 = t748*t744;
  t767 = -1.*t738*t726;
  t769 = t728*t739;
  t770 = t767 + t769;
  t771 = t722*t770;
  t773 = t765 + t771;
  t831 = t726*t728;
  t832 = t738*t739;
  t833 = t831 + t832;
  t834 = t748*t833;
  t839 = t722*t754;
  t840 = t834 + t839;
  t856 = -1.*t722*t754;
  t857 = t765 + t856;
  t867 = -1.*t728*t786;
  t869 = t738*t797;
  t870 = t867 + t869;
  t873 = t738*t786;
  t874 = t728*t797;
  t875 = t873 + t874;
  t885 = t738*t808;
  t886 = -1.*t728*t817;
  t887 = t885 + t886;
  t880 = t728*t808;
  t881 = t738*t817;
  t882 = t880 + t881;
  t859 = -1.*t722*t744;
  t860 = t748*t770;
  t861 = t859 + t860;
  t827 = -1.*t726*t817;
  t829 = t808*t739;
  t830 = t823 + t827 + t828 + t829;
  t853 = -1.*t722*t833;
  t854 = t853 + t755;
  t787 = -1.*t726*t786;
  t799 = t797*t739;
  t801 = t787 + t799;
  t805 = t726*t786;
  t810 = t726*t808;
  t813 = -1.*t797*t739;
  t820 = t817*t739;
  t821 = t805 + t810 + t813 + t820;
  t843 = t823 + t828;
  t855 = 0.14994*t761*t854;
  t858 = 0.14994*t857*t840;
  t863 = 0.14994*t761*t861;
  t864 = 0.14994*t857*t773;
  t865 = t855 + t858 + t863 + t864;
  t871 = -1.*t870*t744;
  t876 = -1.*t875*t754;
  t877 = t871 + t876;
  t879 = t870*t744;
  t883 = t882*t833;
  t884 = t875*t754;
  t888 = t887*t754;
  t889 = t879 + t883 + t884 + t888;
  t891 = t875*t833;
  t892 = t870*t754;
  t893 = t891 + t892;
  t895 = -1.*t875*t744;
  t896 = -1.*t887*t744;
  t897 = -1.*t882*t754;
  t898 = -1.*t870*t770;
  t899 = t895 + t896 + t897 + t898;
  t920 = 0.14994*t857*t877;
  t921 = 0.14994*t857*t889;
  t922 = 0.14994*t893*t861;
  t923 = 0.14994*t854*t899;
  t924 = t920 + t921 + t922 + t923;
  t878 = 0.14994*t761*t877;
  t890 = 0.14994*t761*t889;
  t894 = 0.14994*t893*t773;
  t900 = 0.14994*t840*t899;
  t901 = t878 + t890 + t894 + t900;
  t909 = 0.14994*t830*t854;
  t910 = 0.14994*t801*t857;
  t911 = 0.14994*t821*t857;
  t912 = 0.14994*t843*t861;
  t913 = t909 + t910 + t911 + t912;
  t803 = 0.14994*t801*t761;
  t822 = 0.14994*t821*t761;
  t842 = 0.14994*t830*t840;
  t844 = 0.14994*t843*t773;
  t845 = t803 + t822 + t842 + t844;
  t932 = 0.14994*t821*t893;
  t933 = 0.14994*t830*t877;
  t934 = 0.14994*t843*t889;
  t935 = 0.14994*t801*t899;
  t936 = t932 + t933 + t934 + t935;
  t946 = 0.0033980902199999994*t830;
  t947 = -0.0011052077399999983*t821;
  t948 = t946 + t947;
  t905 = 0.0033980902199999994*t857;
  t906 = -0.0011052077399999983*t861;
  t907 = t905 + t906;
  t762 = 0.0033980902199999994*t761;
  t775 = -0.0011052077399999983*t773;
  t777 = t762 + t775;
  t928 = -0.0011052077399999983*t889;
  t929 = 0.0033980902199999994*t899;
  t930 = t928 + t929;
  p_output1[0]=var2[4]*(-0.5*(0.29988*t761*t773 + 0.29988*t761*t840)*var2[0] - 0.5*t865*var2[1] - 0.5*t901*var2[2] - 0.5*t845*var2[3] - 0.5*t777*var2[4]);
  p_output1[1]=var2[4]*(-0.5*t865*var2[0] - 0.5*(0.29988*t854*t857 + 0.29988*t857*t861)*var2[1] - 0.5*t924*var2[2] - 0.5*t913*var2[3] - 0.5*t907*var2[4]);
  p_output1[2]=var2[4]*(-0.5*t901*var2[0] - 0.5*t924*var2[1] - 0.5*(0.29988*t889*t893 + 0.29988*t877*t899)*var2[2] - 0.5*t936*var2[3] - 0.5*t930*var2[4]);
  p_output1[3]=var2[4]*(-0.5*t845*var2[0] - 0.5*t913*var2[1] - 0.5*t936*var2[2] - 0.5*(0.29988*t801*t830 + 0.29988*t821*t843)*var2[3] - 0.5*t948*var2[4]);
  p_output1[4]=(-0.5*t777*var2[0] - 0.5*t907*var2[1] - 0.5*t930*var2[2] - 0.5*t948*var2[3])*var2[4];
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

#include "Ce2_vec_L4_J5_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L4_J5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
