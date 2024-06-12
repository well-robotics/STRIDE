/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:39:16 GMT-05:00
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
  double t642;
  double t636;
  double t637;
  double t643;
  double t638;
  double t647;
  double t648;
  double t652;
  double t654;
  double t662;
  double t684;
  double t688;
  double t689;
  double t722;
  double t723;
  double t725;
  double t735;
  double t736;
  double t738;
  double t741;
  double t728;
  double t729;
  double t730;
  double t712;
  double t713;
  double t716;
  double t682;
  double t699;
  double t703;
  double t704;
  double t706;
  double t709;
  double t717;
  double t718;
  double t726;
  double t732;
  double t733;
  double t739;
  double t742;
  double t745;
  double t746;
  double t747;
  double t749;
  double t750;
  double t734;
  double t744;
  double t748;
  double t751;
  double t752;
  double t764;
  double t765;
  double t766;
  double t767;
  double t768;
  double t651;
  double t666;
  double t669;
  double t756;
  double t757;
  double t758;
  double t772;
  double t773;
  double t774;
  t642 = Cos(var1[2]);
  t636 = Cos(var1[5]);
  t637 = Sin(var1[2]);
  t643 = Sin(var1[5]);
  t638 = -1.*t636*t637;
  t647 = -1.*t642*t643;
  t648 = t638 + t647;
  t652 = t642*t636;
  t654 = -1.*t637*t643;
  t662 = t652 + t654;
  t684 = t636*t637;
  t688 = t642*t643;
  t689 = t684 + t688;
  t722 = -0.001112*t636;
  t723 = -0.078865*t643;
  t725 = t722 + t723;
  t735 = -0.078865*t636;
  t736 = 0.001112*t643;
  t738 = t735 + t736;
  t741 = t725*t643;
  t728 = 0.001112*t636;
  t729 = 0.078865*t643;
  t730 = t728 + t729;
  t712 = -1.*t642*t636;
  t713 = t637*t643;
  t716 = t712 + t713;
  t682 = 1.28858*t648*t662;
  t699 = Power(t648,2);
  t703 = 0.64429*t699;
  t704 = 0.64429*t648*t689;
  t706 = Power(t662,2);
  t709 = 0.64429*t706;
  t717 = 0.64429*t662*t716;
  t718 = t703 + t704 + t709 + t717;
  t726 = t636*t725;
  t732 = t636*t730;
  t733 = t726 + t732;
  t739 = t636*t738;
  t742 = t739 + t741;
  t745 = -1.*t636*t725;
  t746 = t738*t643;
  t747 = t745 + t746;
  t749 = t730*t643;
  t750 = t741 + t749;
  t734 = 0.64429*t733*t662;
  t744 = 0.64429*t648*t742;
  t748 = 0.64429*t662*t747;
  t751 = 0.64429*t689*t750;
  t752 = t734 + t744 + t748 + t751;
  t764 = 0.64429*t733*t648;
  t765 = 0.64429*t716*t742;
  t766 = 0.64429*t648*t747;
  t767 = 0.64429*t662*t750;
  t768 = t764 + t765 + t766 + t767;
  t651 = -0.050811930850000006*t648;
  t666 = 0.00071645048*t662;
  t669 = t651 + t666;
  t756 = 0.00071645048*t648;
  t757 = -0.050811930850000006*t716;
  t758 = t756 + t757;
  t772 = -0.050811930850000006*t733;
  t773 = 0.00071645048*t750;
  t774 = t772 + t773;
  p_output1[0]=var2[5]*(-0.5*(t682 + 1.28858*t662*t689)*var2[0] - 0.5*t718*var2[1] - 0.5*t752*var2[2] - 0.5*t669*var2[5]);
  p_output1[1]=var2[5]*(-0.5*t718*var2[0] - 0.5*(t682 + 1.28858*t648*t716)*var2[1] - 0.5*t768*var2[2] - 0.5*t758*var2[5]);
  p_output1[2]=var2[5]*(-0.5*t752*var2[0] - 0.5*t768*var2[1] - 0.5*(1.28858*t733*t742 + 1.28858*t747*t750)*var2[2] - 0.5*t774*var2[5]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(-0.5*t669*var2[0] - 0.5*t758*var2[1] - 0.5*t774*var2[2])*var2[5];
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

#include "Ce2_vec_L3_J6_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L3_J6_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
