/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:38:54 GMT-05:00
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
  double t609;
  double t589;
  double t597;
  double t611;
  double t600;
  double t618;
  double t621;
  double t628;
  double t630;
  double t631;
  double t637;
  double t638;
  double t639;
  double t682;
  double t684;
  double t687;
  double t699;
  double t703;
  double t704;
  double t709;
  double t689;
  double t691;
  double t693;
  double t655;
  double t662;
  double t666;
  double t636;
  double t647;
  double t648;
  double t651;
  double t652;
  double t654;
  double t669;
  double t674;
  double t688;
  double t695;
  double t696;
  double t706;
  double t712;
  double t714;
  double t715;
  double t716;
  double t718;
  double t719;
  double t697;
  double t713;
  double t717;
  double t720;
  double t721;
  double t733;
  double t734;
  double t735;
  double t736;
  double t737;
  double t627;
  double t633;
  double t634;
  double t725;
  double t726;
  double t727;
  double t741;
  double t742;
  double t743;
  t609 = Cos(var1[2]);
  t589 = Cos(var1[3]);
  t597 = Sin(var1[2]);
  t611 = Sin(var1[3]);
  t600 = -1.*t589*t597;
  t618 = -1.*t609*t611;
  t621 = t600 + t618;
  t628 = t609*t589;
  t630 = -1.*t597*t611;
  t631 = t628 + t630;
  t637 = t589*t597;
  t638 = t609*t611;
  t639 = t637 + t638;
  t682 = -0.00102*t589;
  t684 = -0.078722*t611;
  t687 = t682 + t684;
  t699 = -0.078722*t589;
  t703 = 0.00102*t611;
  t704 = t699 + t703;
  t709 = t687*t611;
  t689 = 0.00102*t589;
  t691 = 0.078722*t611;
  t693 = t689 + t691;
  t655 = -1.*t609*t589;
  t662 = t597*t611;
  t666 = t655 + t662;
  t636 = 1.29308*t621*t631;
  t647 = Power(t621,2);
  t648 = 0.64654*t647;
  t651 = 0.64654*t621*t639;
  t652 = Power(t631,2);
  t654 = 0.64654*t652;
  t669 = 0.64654*t631*t666;
  t674 = t648 + t651 + t654 + t669;
  t688 = t589*t687;
  t695 = t589*t693;
  t696 = t688 + t695;
  t706 = t589*t704;
  t712 = t706 + t709;
  t714 = -1.*t589*t687;
  t715 = t704*t611;
  t716 = t714 + t715;
  t718 = t693*t611;
  t719 = t709 + t718;
  t697 = 0.64654*t696*t631;
  t713 = 0.64654*t621*t712;
  t717 = 0.64654*t631*t716;
  t720 = 0.64654*t639*t719;
  t721 = t697 + t713 + t717 + t720;
  t733 = 0.64654*t696*t621;
  t734 = 0.64654*t666*t712;
  t735 = 0.64654*t621*t716;
  t736 = 0.64654*t631*t719;
  t737 = t733 + t734 + t735 + t736;
  t627 = -0.05089692188*t621;
  t633 = 0.0006594708000000001*t631;
  t634 = t627 + t633;
  t725 = 0.0006594708000000001*t621;
  t726 = -0.05089692188*t666;
  t727 = t725 + t726;
  t741 = -0.05089692188*t696;
  t742 = 0.0006594708000000001*t719;
  t743 = t741 + t742;
  p_output1[0]=var2[3]*(-0.5*(t636 + 1.29308*t631*t639)*var2[0] - 0.5*t674*var2[1] - 0.5*t721*var2[2] - 0.5*t634*var2[3]);
  p_output1[1]=var2[3]*(-0.5*t674*var2[0] - 0.5*(t636 + 1.29308*t621*t666)*var2[1] - 0.5*t737*var2[2] - 0.5*t727*var2[3]);
  p_output1[2]=var2[3]*(-0.5*t721*var2[0] - 0.5*t737*var2[1] - 0.5*(1.29308*t696*t712 + 1.29308*t716*t719)*var2[2] - 0.5*t743*var2[3]);
  p_output1[3]=(-0.5*t634*var2[0] - 0.5*t727*var2[1] - 0.5*t743*var2[2])*var2[3];
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

#include "Ce2_vec_L2_J4_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L2_J4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
