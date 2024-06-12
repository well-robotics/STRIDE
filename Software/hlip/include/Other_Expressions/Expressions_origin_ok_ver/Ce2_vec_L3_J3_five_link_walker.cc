/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:09 GMT-05:00
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
  double t631;
  double t621;
  double t627;
  double t633;
  double t628;
  double t635;
  double t636;
  double t638;
  double t639;
  double t642;
  double t651;
  double t652;
  double t654;
  double t712;
  double t713;
  double t715;
  double t703;
  double t704;
  double t706;
  double t691;
  double t693;
  double t695;
  double t648;
  double t680;
  double t682;
  double t684;
  double t688;
  double t689;
  double t696;
  double t697;
  double t709;
  double t716;
  double t717;
  double t719;
  double t720;
  double t722;
  double t738;
  double t739;
  double t740;
  double t718;
  double t723;
  double t724;
  double t637;
  double t643;
  double t644;
  double t729;
  double t730;
  double t731;
  t631 = Cos(var1[2]);
  t621 = Cos(var1[5]);
  t627 = Sin(var1[2]);
  t633 = Sin(var1[5]);
  t628 = -1.*t621*t627;
  t635 = -1.*t631*t633;
  t636 = t628 + t635;
  t638 = t631*t621;
  t639 = -1.*t627*t633;
  t642 = t638 + t639;
  t651 = t621*t627;
  t652 = t631*t633;
  t654 = t651 + t652;
  t712 = -0.001112*t621;
  t713 = -0.078865*t633;
  t715 = t712 + t713;
  t703 = -0.078865*t621;
  t704 = 0.001112*t633;
  t706 = t703 + t704;
  t691 = -1.*t631*t621;
  t693 = t627*t633;
  t695 = t691 + t693;
  t648 = 1.28858*t636*t642;
  t680 = Power(t636,2);
  t682 = 0.64429*t680;
  t684 = 0.64429*t636*t654;
  t688 = Power(t642,2);
  t689 = 0.64429*t688;
  t696 = 0.64429*t642*t695;
  t697 = t682 + t684 + t689 + t696;
  t709 = t621*t706;
  t716 = t715*t633;
  t717 = t709 + t716;
  t719 = -1.*t621*t715;
  t720 = t706*t633;
  t722 = t719 + t720;
  t738 = 0.64429*t695*t717;
  t739 = 0.64429*t636*t722;
  t740 = t738 + t739;
  t718 = 0.64429*t636*t717;
  t723 = 0.64429*t642*t722;
  t724 = t718 + t723;
  t637 = -0.050811930850000006*t636;
  t643 = 0.00071645048*t642;
  t644 = t637 + t643;
  t729 = 0.00071645048*t636;
  t730 = -0.050811930850000006*t695;
  t731 = t729 + t730;
  p_output1[0]=var2[2]*(-0.5*(t648 + 1.28858*t642*t654)*var2[0] - 0.5*t697*var2[1] - 0.5*t724*var2[2] - 0.5*t644*var2[5]);
  p_output1[1]=var2[2]*(-0.5*t697*var2[0] - 0.5*(t648 + 1.28858*t636*t695)*var2[1] - 0.5*t740*var2[2] - 0.5*t731*var2[5]);
  p_output1[2]=(-0.5*t724*var2[0] - 0.5*t740*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t644*var2[0] - 0.5*t731*var2[1])*var2[2];
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

#include "Ce2_vec_L3_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L3_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
