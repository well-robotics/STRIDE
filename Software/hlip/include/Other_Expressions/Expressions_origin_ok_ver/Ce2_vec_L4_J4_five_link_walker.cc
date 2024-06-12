/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:28 GMT-05:00
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
  double t703;
  double t682;
  double t693;
  double t704;
  double t726;
  double t699;
  double t716;
  double t720;
  double t662;
  double t722;
  double t728;
  double t729;
  double t735;
  double t736;
  double t738;
  double t759;
  double t760;
  double t761;
  double t762;
  double t765;
  double t766;
  double t754;
  double t755;
  double t756;
  double t741;
  double t742;
  double t744;
  double t746;
  double t747;
  double t748;
  double t781;
  double t783;
  double t785;
  double t786;
  double t787;
  double t790;
  double t800;
  double t801;
  double t812;
  double t813;
  double t814;
  double t816;
  double t817;
  double t818;
  double t824;
  double t825;
  double t826;
  double t804;
  double t805;
  double t806;
  double t757;
  double t767;
  double t769;
  double t771;
  double t773;
  double t775;
  double t796;
  double t797;
  double t799;
  double t803;
  double t807;
  double t808;
  double t809;
  double t815;
  double t819;
  double t820;
  double t822;
  double t823;
  double t827;
  double t828;
  double t829;
  double t831;
  double t832;
  double t833;
  double t835;
  double t836;
  double t837;
  double t838;
  double t839;
  double t858;
  double t859;
  double t860;
  double t861;
  double t862;
  double t821;
  double t830;
  double t834;
  double t840;
  double t841;
  double t849;
  double t850;
  double t851;
  double t770;
  double t777;
  double t778;
  double t870;
  double t871;
  double t872;
  double t845;
  double t846;
  double t847;
  double t739;
  double t749;
  double t751;
  double t866;
  double t867;
  double t868;
  t703 = Cos(var1[3]);
  t682 = Cos(var1[4]);
  t693 = Sin(var1[3]);
  t704 = Sin(var1[4]);
  t726 = Cos(var1[2]);
  t699 = -1.*t682*t693;
  t716 = -1.*t703*t704;
  t720 = t699 + t716;
  t662 = Sin(var1[2]);
  t722 = t662*t720;
  t728 = t703*t682;
  t729 = -1.*t693*t704;
  t735 = t728 + t729;
  t736 = t726*t735;
  t738 = t722 + t736;
  t759 = -1.*t682;
  t760 = 1. + t759;
  t761 = -0.16*t760;
  t762 = -0.167371*t682;
  t765 = 0.022663*t704;
  t766 = t761 + t762 + t765;
  t754 = -0.022663*t682;
  t755 = -0.007370999999999989*t704;
  t756 = t754 + t755;
  t741 = t726*t720;
  t742 = -1.*t703*t682;
  t744 = t693*t704;
  t746 = t742 + t744;
  t747 = t662*t746;
  t748 = t741 + t747;
  t781 = t682*t693;
  t783 = t703*t704;
  t785 = t781 + t783;
  t786 = t726*t785;
  t787 = t662*t735;
  t790 = t786 + t787;
  t800 = -1.*t662*t735;
  t801 = t741 + t800;
  t812 = -1.*t693*t756;
  t813 = t703*t766;
  t814 = t812 + t813;
  t816 = t703*t756;
  t817 = t693*t766;
  t818 = t816 + t817;
  t824 = -1.*t703*t756;
  t825 = -1.*t693*t766;
  t826 = t824 + t825;
  t804 = -1.*t662*t720;
  t805 = t726*t746;
  t806 = t804 + t805;
  t757 = -1.*t682*t756;
  t767 = t766*t704;
  t769 = t757 + t767;
  t771 = t682*t766;
  t773 = t756*t704;
  t775 = t771 + t773;
  t796 = -1.*t662*t785;
  t797 = t796 + t736;
  t799 = 0.14994*t738*t797;
  t803 = 0.14994*t801*t790;
  t807 = 0.14994*t738*t806;
  t808 = 0.14994*t801*t748;
  t809 = t799 + t803 + t807 + t808;
  t815 = -1.*t814*t720;
  t819 = -1.*t818*t735;
  t820 = t815 + t819;
  t822 = t814*t720;
  t823 = t814*t785;
  t827 = t826*t735;
  t828 = t818*t735;
  t829 = t822 + t823 + t827 + t828;
  t831 = t818*t785;
  t832 = t814*t735;
  t833 = t831 + t832;
  t835 = -1.*t826*t720;
  t836 = -1.*t818*t720;
  t837 = -1.*t814*t735;
  t838 = -1.*t814*t746;
  t839 = t835 + t836 + t837 + t838;
  t858 = 0.14994*t801*t820;
  t859 = 0.14994*t801*t829;
  t860 = 0.14994*t833*t806;
  t861 = 0.14994*t797*t839;
  t862 = t858 + t859 + t860 + t861;
  t821 = 0.14994*t738*t820;
  t830 = 0.14994*t738*t829;
  t834 = 0.14994*t833*t748;
  t840 = 0.14994*t790*t839;
  t841 = t821 + t830 + t834 + t840;
  t849 = 0.14994*t769*t801;
  t850 = 0.14994*t775*t806;
  t851 = t849 + t850;
  t770 = 0.14994*t769*t738;
  t777 = 0.14994*t775*t748;
  t778 = t770 + t777;
  t870 = 0.14994*t775*t829;
  t871 = 0.14994*t769*t839;
  t872 = t870 + t871;
  t845 = 0.0033980902199999994*t801;
  t846 = -0.0011052077399999983*t806;
  t847 = t845 + t846;
  t739 = 0.0033980902199999994*t738;
  t749 = -0.0011052077399999983*t748;
  t751 = t739 + t749;
  t866 = -0.0011052077399999983*t829;
  t867 = 0.0033980902199999994*t839;
  t868 = t866 + t867;
  p_output1[0]=var2[3]*(-0.5*(0.29988*t738*t748 + 0.29988*t738*t790)*var2[0] - 0.5*t809*var2[1] - 0.5*t841*var2[2] - 0.5*t778*var2[3] - 0.5*t751*var2[4]);
  p_output1[1]=var2[3]*(-0.5*t809*var2[0] - 0.5*(0.29988*t797*t801 + 0.29988*t801*t806)*var2[1] - 0.5*t862*var2[2] - 0.5*t851*var2[3] - 0.5*t847*var2[4]);
  p_output1[2]=var2[3]*(-0.5*t841*var2[0] - 0.5*t862*var2[1] - 0.5*(0.29988*t829*t833 + 0.29988*t820*t839)*var2[2] - 0.5*t872*var2[3] - 0.5*t868*var2[4]);
  p_output1[3]=(-0.5*t778*var2[0] - 0.5*t851*var2[1] - 0.5*t872*var2[2])*var2[3];
  p_output1[4]=(-0.5*t751*var2[0] - 0.5*t847*var2[1] - 0.5*t868*var2[2])*var2[3];
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

#include "Ce2_vec_L4_J4_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L4_J4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
