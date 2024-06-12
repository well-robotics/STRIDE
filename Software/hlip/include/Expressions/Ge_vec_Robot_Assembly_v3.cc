/*
 * Automatically Generated from Mathematica.
 * Wed 21 Feb 2024 16:40:55 GMT-06:00
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
  double t1594;
  double t1603;
  double t1649;
  double t1671;
  double t1678;
  double t1684;
  double t1683;
  double t1685;
  double t1692;
  double t1698;
  double t1707;
  double t1717;
  double t1722;
  double t1723;
  double t1735;
  double t1736;
  double t1740;
  double t1724;
  double t1727;
  double t1728;
  double t1748;
  double t1749;
  double t1756;
  double t1758;
  double t1815;
  double t1816;
  double t1819;
  double t1821;
  double t1820;
  double t1823;
  double t1825;
  double t1827;
  double t1831;
  double t1834;
  double t1835;
  double t1841;
  double t1842;
  double t1843;
  double t1844;
  double t1847;
  double t1848;
  double t1849;
  double t1854;
  double t1855;
  double t1857;
  double t1859;
  double t1881;
  double t1883;
  double t1885;
  double t1886;
  double t1887;
  double t1889;
  double t1892;
  double t1899;
  double t1901;
  double t1902;
  double t1757;
  double t1761;
  double t1763;
  double t1893;
  double t1895;
  double t1897;
  double t1767;
  double t1768;
  double t1769;
  double t1908;
  double t1913;
  double t1845;
  double t1850;
  double t1934;
  double t1935;
  double t1936;
  double t1937;
  double t1938;
  double t1939;
  double t1940;
  double t1858;
  double t1860;
  double t1861;
  double t1862;
  double t1863;
  double t1864;
  double t1865;
  double t1866;
  double t1867;
  double t1869;
  double t1871;
  double t1872;
  double t1874;
  double t1875;
  double t1876;
  double t1877;
  double t1946;
  double t1947;
  double t1948;
  t1594 = Cos(var1[2]);
  t1603 = Sin(var1[2]);
  t1649 = Cos(var1[3]);
  t1671 = -1.*t1649;
  t1678 = 1. + t1671;
  t1684 = Sin(var1[3]);
  t1683 = 0.017685*t1678;
  t1685 = -0.0695*t1684;
  t1692 = t1683 + t1685;
  t1698 = -1.*t1594*t1692;
  t1707 = -0.0695*t1678;
  t1717 = -0.017685*t1684;
  t1722 = t1707 + t1717;
  t1723 = -1.*t1603*t1722;
  t1735 = -1.*t1594*t1649;
  t1736 = -1.*t1603*t1684;
  t1740 = t1735 + t1736;
  t1724 = -1.*t1649*t1603;
  t1727 = t1594*t1684;
  t1728 = t1724 + t1727;
  t1748 = Cos(var1[4]);
  t1749 = -1.*t1748;
  t1756 = 1. + t1749;
  t1758 = Sin(var1[4]);
  t1815 = Cos(var1[5]);
  t1816 = -1.*t1815;
  t1819 = 1. + t1816;
  t1821 = Sin(var1[5]);
  t1820 = -0.0695*t1819;
  t1823 = 0.017685*t1821;
  t1825 = t1820 + t1823;
  t1827 = -1.*t1603*t1825;
  t1831 = 0.017685*t1819;
  t1834 = 0.0695*t1821;
  t1835 = t1831 + t1834;
  t1841 = -1.*t1594*t1835;
  t1842 = -1.*t1815*t1603;
  t1843 = -1.*t1594*t1821;
  t1844 = t1842 + t1843;
  t1847 = -1.*t1594*t1815;
  t1848 = t1603*t1821;
  t1849 = t1847 + t1848;
  t1854 = Cos(var1[6]);
  t1855 = -1.*t1854;
  t1857 = 1. + t1855;
  t1859 = Sin(var1[6]);
  t1881 = -0.017685*t1649;
  t1883 = t1881 + t1685;
  t1885 = t1594*t1883;
  t1886 = -0.0695*t1649;
  t1887 = 0.017685*t1684;
  t1889 = t1886 + t1887;
  t1892 = -1.*t1603*t1889;
  t1899 = t1594*t1649;
  t1901 = t1603*t1684;
  t1902 = t1899 + t1901;
  t1757 = 0.044273*t1756;
  t1761 = -0.235383*t1758;
  t1763 = t1757 + t1761;
  t1893 = t1649*t1603;
  t1895 = -1.*t1594*t1684;
  t1897 = t1893 + t1895;
  t1767 = -0.235383*t1756;
  t1768 = -0.044273*t1758;
  t1769 = t1767 + t1768;
  t1908 = t1748*t1902;
  t1913 = -1.*t1902*t1758;
  t1845 = -0.15296*t1844;
  t1850 = 0.013861*t1849;
  t1934 = 0.017685*t1815;
  t1935 = -0.0695*t1821;
  t1936 = t1934 + t1935;
  t1937 = t1594*t1936;
  t1938 = 0.0695*t1815;
  t1939 = t1938 + t1823;
  t1940 = -1.*t1603*t1939;
  t1858 = -0.235612*t1857;
  t1860 = 0.042799*t1859;
  t1861 = t1858 + t1860;
  t1862 = t1844*t1861;
  t1863 = 0.042799*t1857;
  t1864 = 0.235612*t1859;
  t1865 = t1863 + t1864;
  t1866 = t1849*t1865;
  t1867 = t1854*t1849;
  t1869 = -1.*t1844*t1859;
  t1871 = t1867 + t1869;
  t1872 = 0.034542*t1871;
  t1874 = t1854*t1844;
  t1875 = t1849*t1859;
  t1876 = t1874 + t1875;
  t1877 = -0.284649*t1876;
  t1946 = t1594*t1815;
  t1947 = -1.*t1603*t1821;
  t1948 = t1946 + t1947;
  p_output1[0]=0;
  p_output1[1]=-30.7727928;
  p_output1[2]=-14.527629000000001*(-0.000041*t1594 + 0.049294*t1603) - 6.5709342*(t1698 + t1723 - 0.152808*t1728 + 0.014762*t1740) - 1.5516477000000002*(t1698 + t1723 + 0.036018*(t1740*t1748 + t1728*t1758) - 0.284422*(t1728*t1748 - 1.*t1740*t1758) + t1740*t1763 + t1728*t1769) - 6.5709342*(t1827 + t1841 + t1845 + t1850) - 1.5516477000000002*(t1827 + t1841 + t1862 + t1866 + t1872 + t1877);
  p_output1[3]=-6.5709342*(t1885 + t1892 - 0.152808*t1897 + 0.014762*t1902) - 1.5516477000000002*(t1885 + t1892 + t1769*t1897 + t1763*t1902 + 0.036018*(t1758*t1897 + t1908) - 0.284422*(t1748*t1897 + t1913));
  p_output1[4]=-1.5516477000000002*(t1728*(-0.235383*t1748 + 0.044273*t1758) + (-0.044273*t1748 + t1761)*t1902 + 0.036018*(-1.*t1728*t1758 + t1908) - 0.284422*(-1.*t1728*t1748 + t1913));
  p_output1[5]=-6.5709342*(t1845 + t1850 + t1937 + t1940) - 1.5516477000000002*(t1862 + t1866 + t1872 + t1877 + t1937 + t1940);
  p_output1[6]=-1.5516477000000002*(t1844*(0.235612*t1854 + t1860) + (0.042799*t1854 - 0.235612*t1859)*t1948 + 0.034542*(t1869 - 1.*t1854*t1948) - 0.284649*(t1874 - 1.*t1859*t1948));
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

#include "Ge_vec_Robot_Assembly_v3.hh"

namespace SymFunction
{

void Ge_vec_Robot_Assembly_v3_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
