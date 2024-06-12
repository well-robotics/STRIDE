/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:38:11 GMT-05:00
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
  double t503;
  double t478;
  double t488;
  double t507;
  double t535;
  double t460;
  double t542;
  double t544;
  double t545;
  double t568;
  double t574;
  double t577;
  double t556;
  double t578;
  double t600;
  double t606;
  double t609;
  double t610;
  double t611;
  double t614;
  double t617;
  double t618;
  double t621;
  double t624;
  double t499;
  double t512;
  double t520;
  double t587;
  double t588;
  double t632;
  double t642;
  double t643;
  double t644;
  double t645;
  double t646;
  double t654;
  double t655;
  double t656;
  double t657;
  double t658;
  double t661;
  double t662;
  double t664;
  double t666;
  double t668;
  double t669;
  double t673;
  double t677;
  double t678;
  double t679;
  double t682;
  double t684;
  double t685;
  double t674;
  double t686;
  double t697;
  double t699;
  double t700;
  double t704;
  double t705;
  double t690;
  double t709;
  double t710;
  double t711;
  double t703;
  double t706;
  double t707;
  double t692;
  t503 = Cos(var1[5]);
  t478 = Cos(var1[6]);
  t488 = Sin(var1[5]);
  t507 = Sin(var1[6]);
  t535 = Cos(var1[2]);
  t460 = Sin(var1[2]);
  t542 = t503*t478;
  t544 = -1.*t488*t507;
  t545 = t542 + t544;
  t568 = -1.*t478*t488;
  t574 = -1.*t503*t507;
  t577 = t568 + t574;
  t556 = t535*t545;
  t578 = t535*t577;
  t600 = t460*t577;
  t606 = t600 + t556;
  t609 = 0.0033974904599999994*t606;
  t610 = -1.*t503*t478;
  t611 = t488*t507;
  t614 = t610 + t611;
  t617 = t460*t614;
  t618 = t578 + t617;
  t621 = -0.0011047579199999977*t618;
  t624 = t609 + t621;
  t499 = t478*t488;
  t512 = t503*t507;
  t520 = t499 + t512;
  t587 = -1.*t460*t545;
  t588 = t578 + t587;
  t632 = -1.*t460*t577;
  t642 = 0.0033974904599999994*t588;
  t643 = t535*t614;
  t644 = t632 + t643;
  t645 = -0.0011047579199999977*t644;
  t646 = t642 + t645;
  t654 = -0.022659*t478;
  t655 = -0.007367999999999986*t507;
  t656 = t654 + t655;
  t657 = -1.*t488*t656;
  t658 = -1.*t478;
  t661 = 1. + t658;
  t662 = -0.16*t661;
  t664 = -0.167368*t478;
  t666 = 0.022659*t507;
  t668 = t662 + t664 + t666;
  t669 = t503*t668;
  t673 = t657 + t669;
  t677 = -1.*t503*t656;
  t678 = -1.*t488*t668;
  t679 = t677 + t678;
  t682 = t503*t656;
  t684 = t488*t668;
  t685 = t682 + t684;
  t674 = t673*t577;
  t686 = t685*t545;
  t697 = 0.022659*t478;
  t699 = 0.007367999999999986*t507;
  t700 = t697 + t699;
  t704 = -0.007367999999999986*t478;
  t705 = t704 + t666;
  t690 = -1.*t685*t577;
  t709 = t503*t700;
  t710 = -1.*t488*t705;
  t711 = t709 + t710;
  t703 = t488*t700;
  t706 = t503*t705;
  t707 = t703 + t706;
  t692 = -1.*t673*t614;
  p_output1[0]=var2[6]*(-0.5*(0.0033974904599999994*(-1.*t460*t520 + t556) - 0.0011047579199999977*t588)*var2[2] - 0.5*t624*var2[5] - 0.5*t624*var2[6]);
  p_output1[1]=var2[6]*(-0.5*(0.0033974904599999994*(-1.*t520*t535 + t587) - 0.0011047579199999977*(-1.*t535*t545 + t632))*var2[2] - 0.5*t646*var2[5] - 0.5*t646*var2[6]);
  p_output1[2]=var2[6]*(-0.5*(-0.0011047579199999977*(t520*t673 + t674 + t545*t679 + t686) + 0.0033974904599999994*(-1.*t545*t673 - 1.*t577*t679 + t690 + t692))*var2[5] - 0.5*(-0.0011047579199999977*(t674 + t686 + t520*t707 + t545*t711) + 0.0033974904599999994*(t690 + t692 - 1.*t545*t707 - 1.*t577*t711))*var2[6]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=-0.5*(0.0033974904599999994*(t507*t656 + t478*t668 + t507*t700 - 1.*t478*t705) - 0.0011047579199999977*(t478*t656 - 1.*t507*t668 + t478*t700 + t507*t705))*Power(var2[6],2);
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

#include "Ce1_vec_L5_J7_five_link_walker.hh"

namespace SymFunction
{

void Ce1_vec_L5_J7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
