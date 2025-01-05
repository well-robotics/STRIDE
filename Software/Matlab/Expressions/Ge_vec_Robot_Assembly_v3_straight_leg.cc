/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:34:49 GMT-05:00
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
  double t1558;
  double t1576;
  double t1624;
  double t1631;
  double t1640;
  double t1647;
  double t1641;
  double t1651;
  double t1661;
  double t1667;
  double t1668;
  double t1669;
  double t1670;
  double t1677;
  double t1689;
  double t1692;
  double t1693;
  double t1678;
  double t1679;
  double t1681;
  double t1705;
  double t1707;
  double t1709;
  double t1721;
  double t1761;
  double t1763;
  double t1764;
  double t1766;
  double t1765;
  double t1768;
  double t1770;
  double t1771;
  double t1772;
  double t1773;
  double t1774;
  double t1775;
  double t1776;
  double t1777;
  double t1778;
  double t1780;
  double t1781;
  double t1782;
  double t1787;
  double t1788;
  double t1789;
  double t1793;
  double t1816;
  double t1817;
  double t1818;
  double t1819;
  double t1820;
  double t1821;
  double t1822;
  double t1827;
  double t1828;
  double t1829;
  double t1711;
  double t1722;
  double t1727;
  double t1823;
  double t1824;
  double t1825;
  double t1732;
  double t1738;
  double t1740;
  double t1835;
  double t1840;
  double t1779;
  double t1783;
  double t1861;
  double t1862;
  double t1863;
  double t1864;
  double t1865;
  double t1866;
  double t1867;
  double t1791;
  double t1795;
  double t1796;
  double t1797;
  double t1798;
  double t1799;
  double t1800;
  double t1801;
  double t1802;
  double t1803;
  double t1805;
  double t1807;
  double t1809;
  double t1810;
  double t1811;
  double t1812;
  double t1873;
  double t1874;
  double t1875;
  t1558 = Cos(var1[2]);
  t1576 = Sin(var1[2]);
  t1624 = Cos(var1[3]);
  t1631 = -1.*t1624;
  t1640 = 1. + t1631;
  t1647 = Sin(var1[3]);
  t1641 = -0.0265*t1640;
  t1651 = -0.0695*t1647;
  t1661 = t1641 + t1651;
  t1667 = -1.*t1558*t1661;
  t1668 = -0.0695*t1640;
  t1669 = 0.0265*t1647;
  t1670 = t1668 + t1669;
  t1677 = -1.*t1576*t1670;
  t1689 = -1.*t1558*t1624;
  t1692 = -1.*t1576*t1647;
  t1693 = t1689 + t1692;
  t1678 = -1.*t1624*t1576;
  t1679 = t1558*t1647;
  t1681 = t1678 + t1679;
  t1705 = Cos(var1[4]);
  t1707 = -1.*t1705;
  t1709 = 1. + t1707;
  t1721 = Sin(var1[4]);
  t1761 = Cos(var1[5]);
  t1763 = -1.*t1761;
  t1764 = 1. + t1763;
  t1766 = Sin(var1[5]);
  t1765 = -0.0695*t1764;
  t1768 = -0.0265*t1766;
  t1770 = t1765 + t1768;
  t1771 = -1.*t1576*t1770;
  t1772 = -0.0265*t1764;
  t1773 = 0.0695*t1766;
  t1774 = t1772 + t1773;
  t1775 = -1.*t1558*t1774;
  t1776 = -1.*t1761*t1576;
  t1777 = -1.*t1558*t1766;
  t1778 = t1776 + t1777;
  t1780 = -1.*t1558*t1761;
  t1781 = t1576*t1766;
  t1782 = t1780 + t1781;
  t1787 = Cos(var1[6]);
  t1788 = -1.*t1787;
  t1789 = 1. + t1788;
  t1793 = Sin(var1[6]);
  t1816 = 0.0265*t1624;
  t1817 = t1816 + t1651;
  t1818 = t1558*t1817;
  t1819 = -0.0695*t1624;
  t1820 = -0.0265*t1647;
  t1821 = t1819 + t1820;
  t1822 = -1.*t1576*t1821;
  t1827 = t1558*t1624;
  t1828 = t1576*t1647;
  t1829 = t1827 + t1828;
  t1711 = -0.0265*t1709;
  t1722 = -0.2375*t1721;
  t1727 = t1711 + t1722;
  t1823 = t1624*t1576;
  t1824 = -1.*t1558*t1647;
  t1825 = t1823 + t1824;
  t1732 = -0.2375*t1709;
  t1738 = 0.0265*t1721;
  t1740 = t1732 + t1738;
  t1835 = t1705*t1829;
  t1840 = -1.*t1829*t1721;
  t1779 = -0.15232*t1778;
  t1783 = -0.025413*t1782;
  t1861 = -0.0265*t1761;
  t1862 = -0.0695*t1766;
  t1863 = t1861 + t1862;
  t1864 = t1558*t1863;
  t1865 = 0.0695*t1761;
  t1866 = t1865 + t1768;
  t1867 = -1.*t1576*t1866;
  t1791 = -0.2375*t1789;
  t1795 = -0.0265*t1793;
  t1796 = t1791 + t1795;
  t1797 = t1778*t1796;
  t1798 = -0.0265*t1789;
  t1799 = 0.2375*t1793;
  t1800 = t1798 + t1799;
  t1801 = t1782*t1800;
  t1802 = t1787*t1782;
  t1803 = -1.*t1778*t1793;
  t1805 = t1802 + t1803;
  t1807 = -0.025226*t1805;
  t1809 = t1787*t1778;
  t1810 = t1782*t1793;
  t1811 = t1809 + t1810;
  t1812 = -0.314506*t1811;
  t1873 = t1558*t1761;
  t1874 = -1.*t1576*t1766;
  t1875 = t1873 + t1874;
  p_output1[0]=0;
  p_output1[1]=-32.734989;
  p_output1[2]=-15.270246*(0.026461*t1558 + 0.046589*t1576) - 6.7739031*(t1667 + t1677 - 0.15232*t1681 - 0.025367*t1693) - 1.9584684*(t1667 + t1677 - 0.025229*(t1693*t1705 + t1681*t1721) - 0.314514*(t1681*t1705 - 1.*t1693*t1721) + t1693*t1727 + t1681*t1740) - 6.7739031*(t1771 + t1775 + t1779 + t1783) - 1.9584684*(t1771 + t1775 + t1797 + t1801 + t1807 + t1812);
  p_output1[3]=-6.7739031*(t1818 + t1822 - 0.15232*t1825 - 0.025367*t1829) - 1.9584684*(t1818 + t1822 + t1740*t1825 + t1727*t1829 - 0.025229*(t1721*t1825 + t1835) - 0.314514*(t1705*t1825 + t1840));
  p_output1[4]=-1.9584684*(t1681*(-0.2375*t1705 - 0.0265*t1721) + (0.0265*t1705 + t1722)*t1829 - 0.025229*(-1.*t1681*t1721 + t1835) - 0.314514*(-1.*t1681*t1705 + t1840));
  p_output1[5]=-6.7739031*(t1779 + t1783 + t1864 + t1867) - 1.9584684*(t1797 + t1801 + t1807 + t1812 + t1864 + t1867);
  p_output1[6]=-1.9584684*(t1778*(0.2375*t1787 + t1795) + (-0.0265*t1787 - 0.2375*t1793)*t1875 - 0.025226*(t1803 - 1.*t1787*t1875) - 0.314506*(t1809 - 1.*t1793*t1875));
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

#include "Ge_vec_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ge_vec_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
