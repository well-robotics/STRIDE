/*
 * Automatically Generated from Mathematica.
 * Wed 21 Feb 2024 16:41:36 GMT-06:00
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
  double t1757;
  double t1767;
  double t1785;
  double t1764;
  double t1769;
  double t1772;
  double t1745;
  double t1802;
  double t1803;
  double t1812;
  double t1727;
  double t1861;
  double t1820;
  double t1881;
  double t1898;
  double t1899;
  double t1784;
  double t1834;
  double t1862;
  double t1863;
  double t1866;
  double t1851;
  double t1867;
  double t1875;
  double t1880;
  double t1901;
  double t1906;
  double t1910;
  double t1915;
  double t1907;
  double t1917;
  double t1919;
  double t1942;
  double t1943;
  double t1936;
  double t1937;
  double t1944;
  double t1945;
  double t1946;
  double t1949;
  double t1951;
  double t1952;
  double t1953;
  double t1954;
  double t1958;
  double t1959;
  double t1960;
  double t1961;
  double t1962;
  double t1938;
  double t1941;
  double t1950;
  double t1955;
  double t1956;
  double t1925;
  double t1926;
  double t1928;
  double t1930;
  double t1931;
  double t1932;
  double t1933;
  double t1934;
  double t1935;
  double t1957;
  double t1963;
  double t1964;
  double t1966;
  double t1967;
  double t1968;
  double t1971;
  double t1972;
  double t1973;
  double t1974;
  double t1975;
  double t1977;
  double t1978;
  double t1979;
  double t1980;
  double t1981;
  double t1965;
  double t1969;
  double t1970;
  double t1992;
  double t1993;
  double t1994;
  double t1976;
  double t1982;
  double t1983;
  double t1995;
  double t1996;
  double t1997;
  double t2006;
  double t2007;
  double t2008;
  double t1984;
  double t1985;
  double t1986;
  double t1998;
  double t1999;
  double t2000;
  double t2009;
  double t2010;
  double t2011;
  double t2017;
  double t2018;
  double t2019;
  t1757 = Cos(var1[6]);
  t1767 = Sin(var1[6]);
  t1785 = Cos(var1[5]);
  t1764 = -0.262015*t1757;
  t1769 = -0.965064*t1767;
  t1772 = t1764 + t1769;
  t1745 = Sin(var1[5]);
  t1802 = 0.965064*t1757;
  t1803 = -0.262015*t1767;
  t1812 = t1802 + t1803;
  t1727 = Cos(var1[2]);
  t1861 = Sin(var1[2]);
  t1820 = t1785*t1812;
  t1881 = 0.262015*t1757;
  t1898 = 0.965064*t1767;
  t1899 = t1881 + t1898;
  t1784 = t1745*t1772;
  t1834 = t1784 + t1820;
  t1862 = t1785*t1772;
  t1863 = -1.*t1745*t1812;
  t1866 = t1862 + t1863;
  t1851 = t1727*t1834;
  t1867 = t1861*t1866;
  t1875 = t1851 + t1867;
  t1880 = t1745*t1812;
  t1901 = t1785*t1899;
  t1906 = t1880 + t1901;
  t1910 = -1.*t1745*t1899;
  t1915 = t1820 + t1910;
  t1907 = t1727*t1906;
  t1917 = t1861*t1915;
  t1919 = t1907 + t1917;
  t1942 = -1.*t1757;
  t1943 = 1. + t1942;
  t1936 = -1.*t1785;
  t1937 = 1. + t1936;
  t1944 = 0.042799*t1943;
  t1945 = 0.034542*t1757;
  t1946 = -0.049037*t1767;
  t1949 = t1944 + t1945 + t1946;
  t1951 = -0.235612*t1943;
  t1952 = -0.284649*t1757;
  t1953 = 0.008256999999999994*t1767;
  t1954 = t1951 + t1952 + t1953;
  t1958 = 0.017685*t1937;
  t1959 = 0.0695*t1745;
  t1960 = t1785*t1949;
  t1961 = t1745*t1954;
  t1962 = t1958 + t1959 + t1960 + t1961;
  t1938 = -0.0695*t1937;
  t1941 = 0.017685*t1745;
  t1950 = -1.*t1745*t1949;
  t1955 = t1785*t1954;
  t1956 = t1938 + t1941 + t1950 + t1955;
  t1925 = -1.*t1861*t1834;
  t1926 = t1727*t1866;
  t1928 = t1925 + t1926;
  t1930 = 0.15817*t1928*t1875;
  t1931 = -1.*t1861*t1906;
  t1932 = t1727*t1915;
  t1933 = t1931 + t1932;
  t1934 = 0.15817*t1933*t1919;
  t1935 = t1930 + t1934;
  t1957 = -1.*t1866*t1956;
  t1963 = -1.*t1834*t1962;
  t1964 = t1957 + t1963;
  t1966 = t1962*t1906;
  t1967 = t1956*t1915;
  t1968 = t1966 + t1967;
  t1971 = 0.017685*t1812;
  t1972 = -1.*t1812*t1949;
  t1973 = -1.*t1772*t1954;
  t1974 = 0.0695*t1899;
  t1975 = t1971 + t1972 + t1973 + t1974;
  t1977 = 0.017685*t1772;
  t1978 = 0.0695*t1812;
  t1979 = t1812*t1954;
  t1980 = t1949*t1899;
  t1981 = t1977 + t1978 + t1979 + t1980;
  t1965 = 0.15817*t1964*t1919;
  t1969 = 0.15817*t1875*t1968;
  t1970 = t1965 + t1969;
  t1992 = 0.15817*t1964*t1933;
  t1993 = 0.15817*t1928*t1968;
  t1994 = t1992 + t1993;
  t1976 = 0.15817*t1919*t1975;
  t1982 = 0.15817*t1875*t1981;
  t1983 = t1976 + t1982;
  t1995 = 0.15817*t1933*t1975;
  t1996 = 0.15817*t1928*t1981;
  t1997 = t1995 + t1996;
  t2006 = 0.15817*t1964*t1975;
  t2007 = 0.15817*t1968*t1981;
  t2008 = 0.000092 + t2006 + t2007;
  t1984 = -0.007827406434441908*t1875;
  t1985 = -0.0007718531672441923*t1919;
  t1986 = t1984 + t1985;
  t1998 = -0.007827406434441908*t1928;
  t1999 = -0.0007718531672441923*t1933;
  t2000 = t1998 + t1999;
  t2009 = -0.0007718531672441923*t1964;
  t2010 = -0.007827406434441908*t1968;
  t2011 = 0.000092 + t2009 + t2010;
  t2017 = -0.0007718531672441923*t1975;
  t2018 = -0.007827406434441908*t1981;
  t2019 = 0.000092 + t2017 + t2018;
  p_output1[0]=0.15817*Power(t1875,2) + 0.15817*Power(t1919,2);
  p_output1[1]=t1935;
  p_output1[2]=t1970;
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=t1983;
  p_output1[6]=t1986;
  p_output1[7]=t1935;
  p_output1[8]=0.15817*Power(t1928,2) + 0.15817*Power(t1933,2);
  p_output1[9]=t1994;
  p_output1[10]=0;
  p_output1[11]=0;
  p_output1[12]=t1997;
  p_output1[13]=t2000;
  p_output1[14]=t1970;
  p_output1[15]=t1994;
  p_output1[16]=0.000092 + 0.15817*Power(t1964,2) + 0.15817*Power(t1968,2);
  p_output1[17]=0;
  p_output1[18]=0;
  p_output1[19]=t2008;
  p_output1[20]=t2011;
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
  p_output1[35]=t1983;
  p_output1[36]=t1997;
  p_output1[37]=t2008;
  p_output1[38]=0;
  p_output1[39]=0;
  p_output1[40]=0.000092 + 0.15817*Power(t1975,2) + 0.15817*Power(t1981,2);
  p_output1[41]=t2019;
  p_output1[42]=t1986;
  p_output1[43]=t2000;
  p_output1[44]=t2011;
  p_output1[45]=0;
  p_output1[46]=0;
  p_output1[47]=t2019;
  p_output1[48]=0.0004831237832820856;
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

#include "Mmat5_Robot_Assembly_v3.hh"

namespace SymFunction
{

void Mmat5_Robot_Assembly_v3_raw(double *p_output1, const double *var1)
{
  // Call Subroutines
  output1(p_output1, var1);

}

}

#endif // MATLAB_MEX_FILE
