/*
 * Automatically Generated from Mathematica.
 * Wed 21 Feb 2024 16:41:35 GMT-06:00
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
  double t1707;
  double t1745;
  double t1683;
  double t1769;
  double t1772;
  double t1775;
  double t1767;
  double t1727;
  double t1757;
  double t1761;
  double t1764;
  double t1784;
  double t1785;
  double t1812;
  double t1820;
  double t1834;
  double t1848;
  double t1851;
  double t1861;
  double t1886;
  double t1887;
  double t1898;
  double t1899;
  double t1901;
  double t1904;
  double t1906;
  double t1907;
  double t1910;
  double t1912;
  double t1867;
  double t1871;
  double t1875;
  double t1876;
  double t1877;
  double t1879;
  double t1880;
  double t1881;
  double t1905;
  double t1914;
  double t1915;
  double t1919;
  double t1920;
  double t1921;
  double t1917;
  double t1924;
  double t1925;
  double t1936;
  double t1937;
  double t1938;
  double t1926;
  double t1928;
  double t1930;
  double t1941;
  double t1942;
  double t1943;
  double t1951;
  double t1952;
  double t1953;
  t1707 = Cos(var1[5]);
  t1745 = Sin(var1[5]);
  t1683 = Sin(var1[2]);
  t1769 = 0.988764*t1707;
  t1772 = 0.149488*t1745;
  t1775 = t1769 + t1772;
  t1767 = Cos(var1[2]);
  t1727 = 0.149488*t1707;
  t1757 = -0.988764*t1745;
  t1761 = t1727 + t1757;
  t1764 = t1683*t1761;
  t1784 = t1767*t1775;
  t1785 = t1764 + t1784;
  t1812 = t1683*t1775;
  t1820 = -0.149488*t1707;
  t1834 = 0.988764*t1745;
  t1848 = t1820 + t1834;
  t1851 = t1767*t1848;
  t1861 = t1812 + t1851;
  t1886 = -1.*t1707;
  t1887 = 1. + t1886;
  t1898 = -0.0695*t1887;
  t1899 = -0.15296*t1707;
  t1901 = 0.0038239999999999993*t1745;
  t1904 = t1898 + t1899 + t1901;
  t1906 = 0.017685*t1887;
  t1907 = 0.013861*t1707;
  t1910 = -0.08346*t1745;
  t1912 = t1906 + t1907 + t1910;
  t1867 = t1767*t1761;
  t1871 = -1.*t1683*t1775;
  t1875 = t1867 + t1871;
  t1876 = 0.66982*t1785*t1875;
  t1877 = -1.*t1683*t1848;
  t1879 = t1784 + t1877;
  t1880 = 0.66982*t1861*t1879;
  t1881 = t1876 + t1880;
  t1905 = -1.*t1761*t1904;
  t1914 = -1.*t1912*t1775;
  t1915 = t1905 + t1914;
  t1919 = t1904*t1775;
  t1920 = t1912*t1848;
  t1921 = t1919 + t1920;
  t1917 = 0.66982*t1915*t1861;
  t1924 = 0.66982*t1785*t1921;
  t1925 = t1917 + t1924;
  t1936 = 0.66982*t1915*t1879;
  t1937 = 0.66982*t1875*t1921;
  t1938 = t1936 + t1937;
  t1926 = -0.05489215178152096*t1785;
  t1928 = 0.010889466036357117*t1861;
  t1930 = t1926 + t1928;
  t1941 = -0.05489215178152096*t1875;
  t1942 = 0.010889466036357117*t1879;
  t1943 = t1941 + t1942;
  t1951 = 0.010889466036357117*t1915;
  t1952 = -0.05489215178152096*t1921;
  t1953 = 0.000248 + t1951 + t1952;
  p_output1[0]=0.66982*Power(t1785,2) + 0.66982*Power(t1861,2);
  p_output1[1]=t1881;
  p_output1[2]=t1925;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t1930;
  p_output1[6]=0;
  p_output1[7]=t1881;
  p_output1[8]=0.66982*Power(t1875,2) + 0.66982*Power(t1879,2);
  p_output1[9]=t1938;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t1943;
  p_output1[13]=0;
  p_output1[14]=t1925;
  p_output1[15]=t1938;
  p_output1[16]=0.000248 + 0.66982*Power(t1915,2) + 0.66982*Power(t1921,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t1953;
  p_output1[20]=0;
  p_output1[21]=0;
  p_output1[22]=0;
  p_output1[23]=0;
  p_output1[24]=0;
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
  p_output1[35]=t1930;
  p_output1[36]=t1943;
  p_output1[37]=t1953;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.004923478184829521;
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

#include "Mmat4_Robot_Assembly_v3.hh"

namespace SymFunction
{

void Mmat4_Robot_Assembly_v3_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
