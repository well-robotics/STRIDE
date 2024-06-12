/*
 * Automatically Generated from Mathematica.
 * Wed 21 Feb 2024 16:41:32 GMT-06:00
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
static void output1(double *p_output1,const double *var1)
{
  double t1623;
  double t1635;
  double t1601;
  double t1727;
  double t1734;
  double t1742;
  double t1717;
  double t1630;
  double t1683;
  double t1685;
  double t1763;
  double t1707;
  double t1745;
  double t1746;
  double t1764;
  double t1767;
  double t1768;
  double t1769;
  double t1772;
  double t1823;
  double t1831;
  double t1851;
  double t1852;
  double t1858;
  double t1860;
  double t1834;
  double t1845;
  double t1848;
  double t1849;
  double t1782;
  double t1784;
  double t1785;
  double t1788;
  double t1802;
  double t1803;
  double t1812;
  double t1820;
  double t1850;
  double t1861;
  double t1862;
  double t1864;
  double t1865;
  double t1866;
  double t1863;
  double t1867;
  double t1871;
  double t1883;
  double t1885;
  double t1886;
  double t1872;
  double t1875;
  double t1876;
  double t1887;
  double t1893;
  double t1897;
  double t1905;
  double t1906;
  double t1907;
  t1623 = Cos(var1[3]);
  t1635 = Sin(var1[3]);
  t1601 = Cos(var1[2]);
  t1727 = 0.987397*t1623;
  t1734 = -0.158262*t1635;
  t1742 = t1727 + t1734;
  t1717 = Sin(var1[2]);
  t1630 = -0.158262*t1623;
  t1683 = -0.987397*t1635;
  t1685 = t1630 + t1683;
  t1763 = t1601*t1742;
  t1707 = t1601*t1685;
  t1745 = t1717*t1742;
  t1746 = t1707 + t1745;
  t1764 = 0.158262*t1623;
  t1767 = 0.987397*t1635;
  t1768 = t1764 + t1767;
  t1769 = t1717*t1768;
  t1772 = t1763 + t1769;
  t1823 = -1.*t1623;
  t1831 = 1. + t1823;
  t1851 = 0.017685*t1831;
  t1852 = 0.014762*t1623;
  t1858 = 0.083308*t1635;
  t1860 = t1851 + t1852 + t1858;
  t1834 = -0.0695*t1831;
  t1845 = -0.152808*t1623;
  t1848 = -0.0029229999999999985*t1635;
  t1849 = t1834 + t1845 + t1848;
  t1782 = -1.*t1717*t1685;
  t1784 = t1782 + t1763;
  t1785 = 0.66982*t1784*t1746;
  t1788 = -1.*t1717*t1742;
  t1802 = t1601*t1768;
  t1803 = t1788 + t1802;
  t1812 = 0.66982*t1803*t1772;
  t1820 = t1785 + t1812;
  t1850 = t1742*t1849;
  t1861 = t1685*t1860;
  t1862 = t1850 + t1861;
  t1864 = -1.*t1742*t1860;
  t1865 = -1.*t1849*t1768;
  t1866 = t1864 + t1865;
  t1863 = 0.66982*t1862*t1772;
  t1867 = 0.66982*t1746*t1866;
  t1871 = t1863 + t1867;
  t1883 = 0.66982*t1862*t1803;
  t1885 = 0.66982*t1784*t1866;
  t1886 = t1883 + t1885;
  t1872 = -0.01076444420770714*t1746;
  t1875 = 0.054788241346999*t1772;
  t1876 = t1872 + t1875;
  t1887 = -0.01076444420770714*t1784;
  t1893 = 0.054788241346999*t1803;
  t1897 = t1887 + t1893;
  t1905 = 0.054788241346999*t1862;
  t1906 = -0.01076444420770714*t1866;
  t1907 = -0.000245 + t1905 + t1906;
  p_output1[0]=0.66982*Power(t1746,2) + 0.66982*Power(t1772,2);
  p_output1[1]=t1820;
  p_output1[2]=t1871;
  p_output1[3]=t1876;
  p_output1[4]=0;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t1820;
  p_output1[8]=0.66982*Power(t1784,2) + 0.66982*Power(t1803,2);
  p_output1[9]=t1886;
  p_output1[10]=t1897;
  p_output1[11]=0;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t1871;
  p_output1[15]=t1886;
  p_output1[16]=0.000245 + 0.66982*Power(t1862,2) + 0.66982*Power(t1866,2);
  p_output1[17]=t1907;
  p_output1[18]=0;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t1876;
  p_output1[22]=t1897;
  p_output1[23]=t1907;
  p_output1[24]=0.004899421559520245;
  p_output1[25]=0;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=0;
  p_output1[29]=0;
  p_output1[30]=0;
  p_output1[31]=0;
  p_output1[32]=0;
  p_output1[33]=0;
  p_output1[34]=0;
  p_output1[35]=0;
  p_output1[36]=0;
  p_output1[37]=0;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0;
  p_output1[41]=0;
  p_output1[42]=0;
  p_output1[43]=0;
  p_output1[44]=0;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=0;
  p_output1[48]=0;
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

  double *var1;
  double *p_output1;

  /*  Check for proper number of arguments.  */ 
  if( nrhs != 1)
    {
      mexErrMsgIdAndTxt("MATLAB:MShaped:invalidNumInputs", "One input(s) required (var1).");
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

  /*  Assign pointers to each input.  */
  var1 = mxGetPr(prhs[0]);
   


   
  /*  Create matrices for return arguments.  */
  plhs[0] = mxCreateDoubleMatrix((mwSize) 7, (mwSize) 7, mxREAL);
  p_output1 = mxGetPr(plhs[0]);


  /* Call the calculation subroutine. */
  output1(p_output1,var1);


}

#else // MATLAB_MEX_FILE

#include "Mmat2_Robot_Assembly_v3.hh"

namespace SymFunction
{

void Mmat2_Robot_Assembly_v3_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
