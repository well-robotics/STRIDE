/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:30:32 GMT-05:00
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
  double t547;
  double t542;
  double t543;
  double t549;
  double t569;
  double t570;
  double t565;
  double t566;
  double t571;
  double t572;
  double t573;
  double t575;
  double t579;
  double t580;
  double t583;
  double t584;
  double t558;
  double t559;
  double t561;
  double t556;
  double t591;
  double t592;
  double t593;
  double t537;
  double t596;
  double t598;
  double t599;
  double t600;
  double t602;
  double t545;
  double t550;
  double t551;
  double t567;
  double t568;
  double t578;
  double t587;
  double t589;
  double t613;
  double t618;
  double t619;
  double t623;
  double t628;
  double t632;
  double t594;
  double t605;
  double t606;
  double t555;
  double t644;
  double t646;
  double t650;
  double t651;
  double t652;
  double t636;
  double t639;
  double t658;
  double t659;
  double t661;
  double t634;
  double t635;
  double t654;
  double t688;
  double t690;
  double t694;
  double t696;
  double t697;
  double t669;
  double t670;
  double t671;
  double t672;
  double t675;
  double t676;
  double t680;
  double t702;
  double t703;
  double t704;
  double t692;
  double t698;
  double t700;
  double t682;
  double t562;
  double t563;
  double t642;
  double t653;
  double t655;
  double t721;
  double t610;
  double t614;
  double t678;
  double t681;
  double t683;
  double t732;
  double t701;
  double t705;
  double t706;
  double t734;
  double t735;
  double t736;
  double t712;
  double t713;
  double t714;
  double t757;
  double t758;
  double t759;
  double t760;
  double t761;
  double t763;
  double t764;
  double t765;
  double t766;
  double t767;
  t547 = Cos(var1[5]);
  t542 = Cos(var1[6]);
  t543 = Sin(var1[5]);
  t549 = Sin(var1[6]);
  t569 = -1.*t542;
  t570 = 1. + t569;
  t565 = -1.*t547;
  t566 = 1. + t565;
  t571 = -0.0265*t570;
  t572 = -0.025226*t542;
  t573 = -0.07700600000000002*t549;
  t575 = t571 + t572 + t573;
  t579 = -0.2375*t570;
  t580 = -0.314506*t542;
  t583 = -0.0012740000000000008*t549;
  t584 = t579 + t580 + t583;
  t558 = t547*t542;
  t559 = -1.*t543*t549;
  t561 = t558 + t559;
  t556 = Sin(var1[2]);
  t591 = t542*t543;
  t592 = t547*t549;
  t593 = t591 + t592;
  t537 = Cos(var1[2]);
  t596 = -0.0695*t566;
  t598 = -0.0265*t543;
  t599 = -1.*t543*t575;
  t600 = t547*t584;
  t602 = t596 + t598 + t599 + t600;
  t545 = -1.*t542*t543;
  t550 = -1.*t547*t549;
  t551 = t545 + t550;
  t567 = -0.0265*t566;
  t568 = 0.0695*t543;
  t578 = t547*t575;
  t587 = t543*t584;
  t589 = t567 + t568 + t578 + t587;
  t613 = t537*t561;
  t618 = -1.*t602*t551;
  t619 = -1.*t589*t561;
  t623 = t618 + t619;
  t628 = t556*t551;
  t632 = t628 + t613;
  t594 = t589*t593;
  t605 = t602*t561;
  t606 = t594 + t605;
  t555 = t537*t551;
  t644 = -0.0265*t547;
  t646 = -0.0695*t543;
  t650 = -1.*t547*t575;
  t651 = -1.*t543*t584;
  t652 = t644 + t646 + t650 + t651;
  t636 = 0.0695*t547;
  t639 = t636 + t598 + t599 + t600;
  t658 = -1.*t547*t542;
  t659 = t543*t549;
  t661 = t658 + t659;
  t634 = 0.19964*t632*t623;
  t635 = t602*t551;
  t654 = t589*t561;
  t688 = -0.07700600000000002*t542;
  t690 = t688 + t583;
  t694 = -0.0012740000000000008*t542;
  t696 = 0.07700600000000002*t549;
  t697 = t694 + t696;
  t669 = t556*t661;
  t670 = t555 + t669;
  t671 = 0.19964*t606*t670;
  t672 = t537*t593;
  t675 = t556*t561;
  t676 = t672 + t675;
  t680 = -1.*t589*t551;
  t702 = -1.*t543*t690;
  t703 = t547*t697;
  t704 = t702 + t703;
  t692 = t547*t690;
  t698 = t543*t697;
  t700 = t692 + t698;
  t682 = -1.*t602*t661;
  t562 = -1.*t556*t561;
  t563 = t555 + t562;
  t642 = t639*t593;
  t653 = t652*t561;
  t655 = t635 + t642 + t653 + t654;
  t721 = -1.*t556*t551;
  t610 = -1.*t556*t593;
  t614 = t610 + t613;
  t678 = -1.*t652*t551;
  t681 = -1.*t639*t561;
  t683 = t678 + t680 + t681 + t682;
  t732 = 0.19964*t563*t623;
  t701 = t700*t593;
  t705 = t704*t561;
  t706 = t635 + t701 + t654 + t705;
  t734 = t537*t661;
  t735 = t721 + t734;
  t736 = 0.19964*t606*t735;
  t712 = -1.*t704*t551;
  t713 = -1.*t700*t561;
  t714 = t680 + t712 + t713 + t682;
  t757 = 0.0695*t542;
  t758 = t542*t584;
  t759 = 0.0265*t549;
  t760 = t575*t549;
  t761 = t757 + t758 + t759 + t760;
  t763 = -0.0265*t542;
  t764 = -1.*t542*t575;
  t765 = 0.0695*t549;
  t766 = t584*t549;
  t767 = t763 + t764 + t765 + t766;
  p_output1[0]=var2[2]*(-0.5*(0.19964*t563*t606 + 0.19964*t614*t623)*var2[2] - 0.5*(t634 + 0.19964*t632*t655 + t671 + 0.19964*t676*t683)*var2[5] - 0.5*(t634 + t671 + 0.19964*t632*t706 + 0.19964*t676*t714)*var2[6]);
  p_output1[1]=var2[2]*(-0.5*(0.19964*(t562 - 1.*t537*t593)*t623 + 0.19964*t606*(-1.*t537*t561 + t721))*var2[2] - 0.5*(0.19964*t563*t655 + 0.19964*t614*t683 + t732 + t736)*var2[5] - 0.5*(0.19964*t563*t706 + 0.19964*t614*t714 + t732 + t736)*var2[6]);
  p_output1[2]=var2[2]*(-0.5*(0.39928*t606*t655 + 0.39928*t623*t683)*var2[5] - 0.5*(0.39928*t606*t706 + 0.39928*t623*t714)*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[2]*(-0.5*(0.19964*t655*t761 + 0.19964*t683*t767)*var2[5] - 0.5*(0.19964*t606*(0.0265*t542 - 0.0695*t549 + t542*t575 - 1.*t549*t584 + t549*t690 + t542*t697) + 0.19964*t623*(-1.*t542*t690 + t549*t697 + t757 + t758 + t759 + t760) + 0.19964*t706*t761 + 0.19964*t714*t767)*var2[6]);
  p_output1[6]=var2[2]*(-0.5*(-0.015373477840000005*t655 - 0.0002543413600000002*t683)*var2[5] - 0.5*(-0.015373477840000005*t706 - 0.0002543413600000002*t714)*var2[6]);
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

#include "Ce1_vec_L5_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce1_vec_L5_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
