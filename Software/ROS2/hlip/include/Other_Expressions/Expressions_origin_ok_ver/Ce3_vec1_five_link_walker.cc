/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:36 GMT-05:00
 */

#ifdef MATLAB_MEX_FILE
#include <stdexcept>
#include <cmath>
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


#endif

#include "mdefs.hpp"
static void output1(double *p_output1,const double *var1,const double *var2)
{
  double t1620;
  double t1617;
  double t1618;
  double t1621;
  double t1669;
  double t1572;
  double t1670;
  double t1676;
  double t1677;
  double t1690;
  double t1696;
  double t1697;
  double t1698;
  double t1699;
  double t1619;
  double t1634;
  double t1667;
  double t1668;
  double t1678;
  double t1679;
  double t1720;
  double t1717;
  double t1718;
  double t1721;
  double t1725;
  double t1726;
  double t1727;
  double t1735;
  double t1736;
  double t1737;
  double t1738;
  double t1739;
  double t1719;
  double t1722;
  double t1723;
  double t1724;
  double t1728;
  double t1729;
  double t1681;
  double t1687;
  double t1688;
  double t1757;
  double t1758;
  double t1759;
  double t1707;
  double t1703;
  double t1704;
  double t1705;
  double t1706;
  double t1708;
  double t1731;
  double t1732;
  double t1733;
  double t1772;
  double t1773;
  double t1774;
  double t1747;
  double t1743;
  double t1744;
  double t1745;
  double t1746;
  double t1748;
  double t1761;
  double t1762;
  double t1763;
  double t1765;
  double t1766;
  double t1768;
  double t1769;
  double t1770;
  double t1776;
  double t1777;
  double t1778;
  double t1780;
  double t1781;
  double t1783;
  double t1784;
  double t1785;
  double t1838;
  double t1839;
  double t1840;
  double t1842;
  double t1843;
  double t1844;
  double t1858;
  double t1859;
  double t1860;
  double t1862;
  double t1863;
  double t1864;
  double t1689;
  double t1700;
  double t1701;
  double t1702;
  double t1710;
  double t1711;
  double t1712;
  double t1713;
  double t1876;
  double t1877;
  double t1878;
  double t1879;
  double t1880;
  double t1760;
  double t1764;
  double t1789;
  double t1790;
  double t1791;
  double t1792;
  double t1793;
  double t1794;
  double t1795;
  double t1796;
  double t1797;
  double t1798;
  double t1832;
  double t1833;
  double t1834;
  double t1835;
  double t1836;
  double t1837;
  double t1841;
  double t1845;
  double t1846;
  double t1848;
  double t1849;
  double t1850;
  double t1899;
  double t1900;
  double t1901;
  double t1881;
  double t1882;
  double t1883;
  double t1886;
  double t1887;
  double t1890;
  double t1891;
  double t1892;
  double t1893;
  double t1894;
  double t1895;
  double t1898;
  double t1903;
  double t1904;
  double t1908;
  double t1933;
  double t1934;
  double t1910;
  double t1936;
  double t1937;
  double t1912;
  double t1734;
  double t1740;
  double t1741;
  double t1742;
  double t1750;
  double t1751;
  double t1752;
  double t1753;
  double t1949;
  double t1950;
  double t1951;
  double t1952;
  double t1953;
  double t1775;
  double t1779;
  double t1809;
  double t1810;
  double t1811;
  double t1812;
  double t1813;
  double t1814;
  double t1815;
  double t1816;
  double t1817;
  double t1818;
  double t1852;
  double t1853;
  double t1854;
  double t1855;
  double t1856;
  double t1857;
  double t1861;
  double t1865;
  double t1866;
  double t1868;
  double t1869;
  double t1870;
  double t1972;
  double t1973;
  double t1974;
  double t1954;
  double t1955;
  double t1956;
  double t1959;
  double t1960;
  double t1963;
  double t1964;
  double t1965;
  double t1966;
  double t1967;
  double t1968;
  double t1971;
  double t1976;
  double t1977;
  double t1981;
  double t2006;
  double t2007;
  double t1983;
  double t2009;
  double t2010;
  double t1985;
  t1620 = Cos(var1[3]);
  t1617 = Cos(var1[4]);
  t1618 = Sin(var1[3]);
  t1621 = Sin(var1[4]);
  t1669 = Sin(var1[2]);
  t1572 = Cos(var1[2]);
  t1670 = t1620*t1617;
  t1676 = -1.*t1618*t1621;
  t1677 = t1670 + t1676;
  t1690 = -1.*t1617;
  t1696 = 1. + t1690;
  t1697 = 0.5*t1696;
  t1698 = 0.671885*t1617;
  t1699 = t1697 + t1698;
  t1619 = -1.*t1617*t1618;
  t1634 = -1.*t1620*t1621;
  t1667 = t1619 + t1634;
  t1668 = t1572*t1667;
  t1678 = -1.*t1669*t1677;
  t1679 = t1668 + t1678;
  t1720 = Cos(var1[5]);
  t1717 = Cos(var1[6]);
  t1718 = Sin(var1[5]);
  t1721 = Sin(var1[6]);
  t1725 = t1720*t1717;
  t1726 = -1.*t1718*t1721;
  t1727 = t1725 + t1726;
  t1735 = -1.*t1717;
  t1736 = 1. + t1735;
  t1737 = 0.5*t1736;
  t1738 = 0.671885*t1717;
  t1739 = t1737 + t1738;
  t1719 = -1.*t1717*t1718;
  t1722 = -1.*t1720*t1721;
  t1723 = t1719 + t1722;
  t1724 = t1572*t1723;
  t1728 = -1.*t1669*t1727;
  t1729 = t1724 + t1728;
  t1681 = -1.*t1620*t1669;
  t1687 = -1.*t1572*t1618;
  t1688 = t1681 + t1687;
  t1757 = t1572*t1620;
  t1758 = -1.*t1669*t1618;
  t1759 = t1757 + t1758;
  t1707 = t1572*t1677;
  t1703 = t1617*t1618;
  t1704 = t1620*t1621;
  t1705 = t1703 + t1704;
  t1706 = -1.*t1669*t1705;
  t1708 = t1706 + t1707;
  t1731 = -1.*t1720*t1669;
  t1732 = -1.*t1572*t1718;
  t1733 = t1731 + t1732;
  t1772 = t1572*t1720;
  t1773 = -1.*t1669*t1718;
  t1774 = t1772 + t1773;
  t1747 = t1572*t1727;
  t1743 = t1717*t1718;
  t1744 = t1720*t1721;
  t1745 = t1743 + t1744;
  t1746 = -1.*t1669*t1745;
  t1748 = t1746 + t1747;
  t1761 = t1620*t1669;
  t1762 = t1572*t1618;
  t1763 = t1761 + t1762;
  t1765 = t1669*t1667;
  t1766 = t1765 + t1707;
  t1768 = t1572*t1705;
  t1769 = t1669*t1677;
  t1770 = t1768 + t1769;
  t1776 = t1720*t1669;
  t1777 = t1572*t1718;
  t1778 = t1776 + t1777;
  t1780 = t1669*t1723;
  t1781 = t1780 + t1747;
  t1783 = t1572*t1745;
  t1784 = t1669*t1727;
  t1785 = t1783 + t1784;
  t1838 = t1699*t1618;
  t1839 = 0.171885*t1620*t1621;
  t1840 = t1838 + t1839;
  t1842 = t1620*t1699;
  t1843 = -0.171885*t1618*t1621;
  t1844 = t1842 + t1843;
  t1858 = t1739*t1718;
  t1859 = 0.171885*t1720*t1721;
  t1860 = t1858 + t1859;
  t1862 = t1720*t1739;
  t1863 = -0.171885*t1718*t1721;
  t1864 = t1862 + t1863;
  t1689 = -0.51185934*t1688;
  t1700 = t1699*t1621;
  t1701 = -0.171885*t1617*t1621;
  t1702 = t1700 + t1701;
  t1710 = t1699*t1617;
  t1711 = Power(t1621,2);
  t1712 = 0.171885*t1711;
  t1713 = t1710 + t1712;
  t1876 = -1.*t1620*t1617;
  t1877 = t1618*t1621;
  t1878 = t1876 + t1877;
  t1879 = t1669*t1878;
  t1880 = t1668 + t1879;
  t1760 = -6.8522*t1688*t1759;
  t1764 = -6.8522*t1763*t1759;
  t1789 = Power(t1688,2);
  t1790 = -3.4261*t1789;
  t1791 = -3.4261*t1688*t1763;
  t1792 = Power(t1759,2);
  t1793 = -3.4261*t1792;
  t1794 = -1.*t1572*t1620;
  t1795 = t1669*t1618;
  t1796 = t1794 + t1795;
  t1797 = -3.4261*t1759*t1796;
  t1798 = -1.*t1669*t1667;
  t1832 = Power(t1620,2);
  t1833 = 0.1494*t1832;
  t1834 = Power(t1618,2);
  t1835 = 0.1494*t1834;
  t1836 = t1833 + t1835;
  t1837 = -3.4261*t1688*t1836;
  t1841 = -1.*t1840*t1677;
  t1845 = -1.*t1667*t1844;
  t1846 = t1841 + t1845;
  t1848 = t1840*t1705;
  t1849 = t1677*t1844;
  t1850 = t1848 + t1849;
  t1899 = -1.*t1699*t1618;
  t1900 = -0.171885*t1620*t1621;
  t1901 = t1899 + t1900;
  t1881 = 0.0732367608*var2[4]*t1880;
  t1882 = -0.85216*t1702*t1766;
  t1883 = -0.85216*t1713*t1880;
  t1886 = -1.70432*t1766*t1770;
  t1887 = -1.70432*t1766*t1880;
  t1890 = -0.85216*t1766*t1708;
  t1891 = -0.85216*t1679*t1770;
  t1892 = t1572*t1878;
  t1893 = t1798 + t1892;
  t1894 = -0.85216*t1766*t1893;
  t1895 = -0.85216*t1679*t1880;
  t1898 = -0.85216*t1766*t1846;
  t1903 = t1840*t1677;
  t1904 = t1667*t1844;
  t1908 = -0.85216*t1850*t1880;
  t1933 = -0.171885*t1617*t1618;
  t1934 = t1933 + t1900;
  t1910 = -1.*t1667*t1840;
  t1936 = 0.171885*t1620*t1617;
  t1937 = t1936 + t1843;
  t1912 = -1.*t1844*t1878;
  t1734 = -0.51185934*t1733;
  t1740 = t1739*t1721;
  t1741 = -0.171885*t1717*t1721;
  t1742 = t1740 + t1741;
  t1750 = t1739*t1717;
  t1751 = Power(t1721,2);
  t1752 = 0.171885*t1751;
  t1753 = t1750 + t1752;
  t1949 = -1.*t1720*t1717;
  t1950 = t1718*t1721;
  t1951 = t1949 + t1950;
  t1952 = t1669*t1951;
  t1953 = t1724 + t1952;
  t1775 = -6.8522*t1733*t1774;
  t1779 = -6.8522*t1778*t1774;
  t1809 = Power(t1733,2);
  t1810 = -3.4261*t1809;
  t1811 = -3.4261*t1733*t1778;
  t1812 = Power(t1774,2);
  t1813 = -3.4261*t1812;
  t1814 = -1.*t1572*t1720;
  t1815 = t1669*t1718;
  t1816 = t1814 + t1815;
  t1817 = -3.4261*t1774*t1816;
  t1818 = -1.*t1669*t1723;
  t1852 = Power(t1720,2);
  t1853 = 0.1494*t1852;
  t1854 = Power(t1718,2);
  t1855 = 0.1494*t1854;
  t1856 = t1853 + t1855;
  t1857 = -3.4261*t1733*t1856;
  t1861 = -1.*t1860*t1727;
  t1865 = -1.*t1723*t1864;
  t1866 = t1861 + t1865;
  t1868 = t1860*t1745;
  t1869 = t1727*t1864;
  t1870 = t1868 + t1869;
  t1972 = -1.*t1739*t1718;
  t1973 = -0.171885*t1720*t1721;
  t1974 = t1972 + t1973;
  t1954 = 0.0732367608*var2[6]*t1953;
  t1955 = -0.85216*t1742*t1781;
  t1956 = -0.85216*t1753*t1953;
  t1959 = -1.70432*t1781*t1785;
  t1960 = -1.70432*t1781*t1953;
  t1963 = -0.85216*t1781*t1748;
  t1964 = -0.85216*t1729*t1785;
  t1965 = t1572*t1951;
  t1966 = t1818 + t1965;
  t1967 = -0.85216*t1781*t1966;
  t1968 = -0.85216*t1729*t1953;
  t1971 = -0.85216*t1781*t1866;
  t1976 = t1860*t1727;
  t1977 = t1723*t1864;
  t1981 = -0.85216*t1870*t1953;
  t2006 = -0.171885*t1717*t1718;
  t2007 = t2006 + t1973;
  t1983 = -1.*t1723*t1860;
  t2009 = 0.171885*t1720*t1717;
  t2010 = t2009 + t1863;
  t1985 = -1.*t1864*t1951;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[0]*(-0.5*(t1760 + t1764 - 1.70432*t1679*t1766 - 1.70432*t1708*t1770 + t1775 + t1779 - 1.70432*t1729*t1781 - 1.70432*t1748*t1785)*var2[0] - 0.5*(-0.85216*Power(t1679,2) - 0.85216*Power(t1708,2) - 0.85216*Power(t1729,2) - 0.85216*Power(t1748,2) - 0.85216*(t1678 - 1.*t1572*t1705)*t1770 - 0.85216*(t1728 - 1.*t1572*t1745)*t1785 + t1790 + t1791 + t1793 + t1797 - 0.85216*t1766*(-1.*t1572*t1677 + t1798) + t1810 + t1811 + t1813 + t1817 - 0.85216*t1781*(-1.*t1572*t1727 + t1818))*var2[1] - 0.5*(3.70591*t1669 + t1837 - 0.85216*t1708*t1846 - 0.85216*t1679*t1850 + t1857 - 0.85216*t1748*t1866 - 0.85216*t1729*t1870)*var2[2] - 0.5*(t1689 - 0.85216*t1702*t1708 - 0.85216*t1679*t1713)*var2[3] + 0.0732367608*t1679*var2[4] - 0.5*(t1734 - 0.85216*t1742*t1748 - 0.85216*t1729*t1753)*var2[5] + 0.0732367608*t1729*var2[6]);
  p_output1[3]=var2[0]*(t1881 - 0.5*(t1760 + t1764 + t1886 + t1887)*var2[0] - 0.5*(t1790 + t1791 + t1793 + t1797 + t1890 + t1891 + t1894 + t1895)*var2[1] - 0.5*(t1837 + t1898 - 0.85216*t1766*(t1705*t1844 + t1677*t1901 + t1903 + t1904) + t1908 - 0.85216*t1770*(-1.*t1677*t1844 - 1.*t1667*t1901 + t1910 + t1912))*var2[2] - 0.5*(t1689 + t1882 + t1883)*var2[3]);
  p_output1[4]=var2[0]*(t1881 - 0.5*(t1886 + t1887)*var2[0] - 0.5*(t1890 + t1891 + t1894 + t1895)*var2[1] - 0.5*(t1898 + t1908 - 0.85216*t1770*(t1910 + t1912 - 1.*t1667*t1934 - 1.*t1677*t1937) - 0.85216*t1766*(t1903 + t1904 + t1677*t1934 + t1705*t1937))*var2[2] - 0.5*(-0.85216*(0.171885*t1617*t1621 - 1.*t1621*t1699)*t1766 - 0.85216*(-0.171885*Power(t1617,2) + t1710)*t1770 + t1882 + t1883)*var2[3]);
  p_output1[5]=var2[0]*(t1954 - 0.5*(t1775 + t1779 + t1959 + t1960)*var2[0] - 0.5*(t1810 + t1811 + t1813 + t1817 + t1963 + t1964 + t1967 + t1968)*var2[1] - 0.5*(t1857 + t1971 - 0.85216*t1781*(t1745*t1864 + t1727*t1974 + t1976 + t1977) + t1981 - 0.85216*t1785*(-1.*t1727*t1864 - 1.*t1723*t1974 + t1983 + t1985))*var2[2] - 0.5*(t1734 + t1955 + t1956)*var2[5]);
  p_output1[6]=var2[0]*(t1954 - 0.5*(t1959 + t1960)*var2[0] - 0.5*(t1963 + t1964 + t1967 + t1968)*var2[1] - 0.5*(t1971 + t1981 - 0.85216*t1785*(t1983 + t1985 - 1.*t1723*t2007 - 1.*t1727*t2010) - 0.85216*t1781*(t1976 + t1977 + t1727*t2007 + t1745*t2010))*var2[2] - 0.5*(-0.85216*(0.171885*t1717*t1721 - 1.*t1721*t1739)*t1781 - 0.85216*(-0.171885*Power(t1717,2) + t1750)*t1785 + t1955 + t1956)*var2[5]);
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

#include "Ce3_vec1_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec1_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
