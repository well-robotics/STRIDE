/*
 * Automatically Generated from Mathematica.
 * Wed 21 Feb 2024 16:41:33 GMT-06:00
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
  double t1635;
  double t1707;
  double t1746;
  double t1683;
  double t1727;
  double t1734;
  double t1630;
  double t1757;
  double t1761;
  double t1763;
  double t1601;
  double t1772;
  double t1764;
  double t1820;
  double t1834;
  double t1845;
  double t1745;
  double t1767;
  double t1769;
  double t1775;
  double t1777;
  double t1781;
  double t1784;
  double t1785;
  double t1812;
  double t1848;
  double t1850;
  double t1851;
  double t1852;
  double t1858;
  double t1861;
  double t1862;
  double t1893;
  double t1897;
  double t1883;
  double t1885;
  double t1898;
  double t1899;
  double t1901;
  double t1903;
  double t1905;
  double t1906;
  double t1907;
  double t1909;
  double t1886;
  double t1887;
  double t1904;
  double t1910;
  double t1911;
  double t1914;
  double t1915;
  double t1916;
  double t1917;
  double t1918;
  double t1871;
  double t1872;
  double t1875;
  double t1876;
  double t1877;
  double t1878;
  double t1879;
  double t1880;
  double t1881;
  double t1912;
  double t1919;
  double t1920;
  double t1922;
  double t1923;
  double t1924;
  double t1927;
  double t1928;
  double t1929;
  double t1930;
  double t1931;
  double t1933;
  double t1934;
  double t1935;
  double t1936;
  double t1937;
  double t1921;
  double t1925;
  double t1926;
  double t1952;
  double t1953;
  double t1954;
  double t1932;
  double t1938;
  double t1941;
  double t1955;
  double t1956;
  double t1957;
  double t1966;
  double t1967;
  double t1968;
  double t1942;
  double t1943;
  double t1944;
  double t1958;
  double t1959;
  double t1960;
  double t1969;
  double t1970;
  double t1971;
  double t1977;
  double t1978;
  double t1979;
  t1635 = Cos(var1[4]);
  t1707 = Sin(var1[4]);
  t1746 = Cos(var1[3]);
  t1683 = 0.262025*t1635;
  t1727 = -0.965061*t1707;
  t1734 = t1683 + t1727;
  t1630 = Sin(var1[3]);
  t1757 = 0.965061*t1635;
  t1761 = 0.262025*t1707;
  t1763 = t1757 + t1761;
  t1601 = Sin(var1[2]);
  t1772 = Cos(var1[2]);
  t1764 = t1746*t1763;
  t1820 = -0.262025*t1635;
  t1834 = 0.965061*t1707;
  t1845 = t1820 + t1834;
  t1745 = t1630*t1734;
  t1767 = t1745 + t1764;
  t1769 = t1601*t1767;
  t1775 = t1746*t1734;
  t1777 = -1.*t1630*t1763;
  t1781 = t1775 + t1777;
  t1784 = t1772*t1781;
  t1785 = t1769 + t1784;
  t1812 = t1630*t1763;
  t1848 = t1746*t1845;
  t1850 = t1812 + t1848;
  t1851 = t1601*t1850;
  t1852 = -1.*t1630*t1845;
  t1858 = t1764 + t1852;
  t1861 = t1772*t1858;
  t1862 = t1851 + t1861;
  t1893 = -1.*t1635;
  t1897 = 1. + t1893;
  t1883 = -1.*t1746;
  t1885 = 1. + t1883;
  t1898 = -0.235383*t1897;
  t1899 = -0.284422*t1635;
  t1901 = -0.008254999999999998*t1707;
  t1903 = t1898 + t1899 + t1901;
  t1905 = 0.044273*t1897;
  t1906 = 0.036018*t1635;
  t1907 = 0.049039*t1707;
  t1909 = t1905 + t1906 + t1907;
  t1886 = -0.0695*t1885;
  t1887 = -0.017685*t1630;
  t1904 = t1746*t1903;
  t1910 = t1630*t1909;
  t1911 = t1886 + t1887 + t1904 + t1910;
  t1914 = 0.017685*t1885;
  t1915 = -0.0695*t1630;
  t1916 = -1.*t1630*t1903;
  t1917 = t1746*t1909;
  t1918 = t1914 + t1915 + t1916 + t1917;
  t1871 = t1772*t1767;
  t1872 = -1.*t1601*t1781;
  t1875 = t1871 + t1872;
  t1876 = 0.15817*t1785*t1875;
  t1877 = t1772*t1850;
  t1878 = -1.*t1601*t1858;
  t1879 = t1877 + t1878;
  t1880 = 0.15817*t1862*t1879;
  t1881 = t1876 + t1880;
  t1912 = t1911*t1767;
  t1919 = t1918*t1781;
  t1920 = t1912 + t1919;
  t1922 = -1.*t1911*t1850;
  t1923 = -1.*t1918*t1858;
  t1924 = t1922 + t1923;
  t1927 = -1.*t1734*t1909;
  t1928 = -0.0695*t1763;
  t1929 = -1.*t1903*t1763;
  t1930 = -0.017685*t1845;
  t1931 = t1927 + t1928 + t1929 + t1930;
  t1933 = -0.0695*t1734;
  t1934 = -0.017685*t1763;
  t1935 = t1909*t1763;
  t1936 = t1903*t1845;
  t1937 = t1933 + t1934 + t1935 + t1936;
  t1921 = 0.15817*t1920*t1862;
  t1925 = 0.15817*t1785*t1924;
  t1926 = t1921 + t1925;
  t1952 = 0.15817*t1920*t1879;
  t1953 = 0.15817*t1875*t1924;
  t1954 = t1952 + t1953;
  t1932 = 0.15817*t1862*t1931;
  t1938 = 0.15817*t1785*t1937;
  t1941 = t1932 + t1938;
  t1955 = 0.15817*t1879*t1931;
  t1956 = 0.15817*t1875*t1937;
  t1957 = t1955 + t1956;
  t1966 = 0.15817*t1920*t1931;
  t1967 = 0.15817*t1924*t1937;
  t1968 = -0.000092 + t1966 + t1967;
  t1942 = 0.0007723228234814011*t1785;
  t1943 = 0.007827618624400175*t1862;
  t1944 = t1942 + t1943;
  t1958 = 0.0007723228234814011*t1875;
  t1959 = 0.007827618624400175*t1879;
  t1960 = t1958 + t1959;
  t1969 = 0.007827618624400175*t1920;
  t1970 = 0.0007723228234814011*t1924;
  t1971 = -0.000092 + t1969 + t1970;
  t1977 = 0.007827618624400175*t1931;
  t1978 = 0.0007723228234814011*t1937;
  t1979 = 0.000092 + t1977 + t1978;
  p_output1[0]=0.15817*Power(t1785,2) + 0.15817*Power(t1862,2);
  p_output1[1]=t1881;
  p_output1[2]=t1926;
  p_output1[3]=t1941;
  p_output1[4]=t1944;
  p_output1[5]=0;
  p_output1[6]=0;
  p_output1[7]=t1881;
  p_output1[8]=0.15817*Power(t1875,2) + 0.15817*Power(t1879,2);
  p_output1[9]=t1954;
  p_output1[10]=t1957;
  p_output1[11]=t1960;
  p_output1[12]=0;
  p_output1[13]=0;
  p_output1[14]=t1926;
  p_output1[15]=t1954;
  p_output1[16]=0.000092 + 0.15817*Power(t1920,2) + 0.15817*Power(t1924,2);
  p_output1[17]=t1968;
  p_output1[18]=t1971;
  p_output1[19]=0;
  p_output1[20]=0;
  p_output1[21]=t1941;
  p_output1[22]=t1957;
  p_output1[23]=t1968;
  p_output1[24]=0.000092 + 0.15817*Power(t1931,2) + 0.15817*Power(t1937,2);
  p_output1[25]=t1979;
  p_output1[26]=0;
  p_output1[27]=0;
  p_output1[28]=t1944;
  p_output1[29]=t1960;
  p_output1[30]=t1971;
  p_output1[31]=t1979;
  p_output1[32]=0.0004831493701253511;
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

#include "Mmat3_Robot_Assembly_v3.hh"

namespace SymFunction
{

void Mmat3_Robot_Assembly_v3_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
