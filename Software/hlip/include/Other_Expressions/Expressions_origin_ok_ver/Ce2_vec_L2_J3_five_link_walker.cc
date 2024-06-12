/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:38:52 GMT-05:00
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
  double t557;
  double t512;
  double t525;
  double t578;
  double t542;
  double t588;
  double t589;
  double t600;
  double t606;
  double t609;
  double t627;
  double t628;
  double t630;
  double t655;
  double t662;
  double t664;
  double t648;
  double t651;
  double t652;
  double t641;
  double t642;
  double t643;
  double t621;
  double t635;
  double t636;
  double t637;
  double t638;
  double t639;
  double t644;
  double t645;
  double t654;
  double t666;
  double t669;
  double t675;
  double t680;
  double t682;
  double t704;
  double t706;
  double t708;
  double t674;
  double t684;
  double t686;
  double t597;
  double t611;
  double t617;
  double t691;
  double t693;
  double t694;
  t557 = Cos(var1[2]);
  t512 = Cos(var1[3]);
  t525 = Sin(var1[2]);
  t578 = Sin(var1[3]);
  t542 = -1.*t512*t525;
  t588 = -1.*t557*t578;
  t589 = t542 + t588;
  t600 = t557*t512;
  t606 = -1.*t525*t578;
  t609 = t600 + t606;
  t627 = t512*t525;
  t628 = t557*t578;
  t630 = t627 + t628;
  t655 = -0.00102*t512;
  t662 = -0.078722*t578;
  t664 = t655 + t662;
  t648 = -0.078722*t512;
  t651 = 0.00102*t578;
  t652 = t648 + t651;
  t641 = -1.*t557*t512;
  t642 = t525*t578;
  t643 = t641 + t642;
  t621 = 1.29308*t589*t609;
  t635 = Power(t589,2);
  t636 = 0.64654*t635;
  t637 = 0.64654*t589*t630;
  t638 = Power(t609,2);
  t639 = 0.64654*t638;
  t644 = 0.64654*t609*t643;
  t645 = t636 + t637 + t639 + t644;
  t654 = t512*t652;
  t666 = t664*t578;
  t669 = t654 + t666;
  t675 = -1.*t512*t664;
  t680 = t652*t578;
  t682 = t675 + t680;
  t704 = 0.64654*t643*t669;
  t706 = 0.64654*t589*t682;
  t708 = t704 + t706;
  t674 = 0.64654*t589*t669;
  t684 = 0.64654*t609*t682;
  t686 = t674 + t684;
  t597 = -0.05089692188*t589;
  t611 = 0.0006594708000000001*t609;
  t617 = t597 + t611;
  t691 = 0.0006594708000000001*t589;
  t693 = -0.05089692188*t643;
  t694 = t691 + t693;
  p_output1[0]=var2[2]*(-0.5*(t621 + 1.29308*t609*t630)*var2[0] - 0.5*t645*var2[1] - 0.5*t686*var2[2] - 0.5*t617*var2[3]);
  p_output1[1]=var2[2]*(-0.5*t645*var2[0] - 0.5*(t621 + 1.29308*t589*t643)*var2[1] - 0.5*t708*var2[2] - 0.5*t694*var2[3]);
  p_output1[2]=(-0.5*t686*var2[0] - 0.5*t708*var2[1])*var2[2];
  p_output1[3]=(-0.5*t617*var2[0] - 0.5*t694*var2[1])*var2[2];
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

#include "Ce2_vec_L2_J3_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec_L2_J3_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
