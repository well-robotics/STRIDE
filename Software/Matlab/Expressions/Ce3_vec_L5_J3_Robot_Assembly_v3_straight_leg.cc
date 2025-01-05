/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:34:07 GMT-05:00
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
  double t1519;
  double t1514;
  double t1517;
  double t1524;
  double t1552;
  double t1553;
  double t1547;
  double t1548;
  double t1554;
  double t1557;
  double t1558;
  double t1559;
  double t1562;
  double t1564;
  double t1566;
  double t1568;
  double t1536;
  double t1537;
  double t1541;
  double t1535;
  double t1573;
  double t1574;
  double t1575;
  double t1508;
  double t1577;
  double t1578;
  double t1579;
  double t1581;
  double t1583;
  double t1518;
  double t1529;
  double t1530;
  double t1549;
  double t1551;
  double t1560;
  double t1569;
  double t1570;
  double t1576;
  double t1584;
  double t1585;
  double t1543;
  double t1598;
  double t1599;
  double t1601;
  double t1640;
  double t1641;
  double t1642;
  double t1645;
  double t1646;
  double t1631;
  double t1632;
  double t1630;
  double t1639;
  double t1647;
  double t1648;
  double t1649;
  double t1656;
  double t1657;
  double t1658;
  double t1660;
  double t1661;
  double t1662;
  double t1666;
  double t1667;
  double t1531;
  double t1544;
  double t1610;
  double t1591;
  double t1596;
  double t1597;
  double t1697;
  double t1698;
  double t1719;
  double t1720;
  double t1722;
  double t1723;
  double t1724;
  double t1728;
  double t1729;
  double t1730;
  double t1721;
  double t1725;
  double t1726;
  double t1671;
  double t1673;
  double t1675;
  double t1676;
  double t1677;
  double t1727;
  double t1731;
  double t1732;
  double t1679;
  double t1680;
  double t1681;
  double t1683;
  double t1684;
  double t1734;
  double t1735;
  double t1736;
  double t1689;
  double t1691;
  double t1692;
  double t1693;
  double t1699;
  double t1701;
  double t1703;
  double t1705;
  double t1707;
  double t1708;
  double t1709;
  t1519 = Cos(var1[5]);
  t1514 = Cos(var1[6]);
  t1517 = Sin(var1[5]);
  t1524 = Sin(var1[6]);
  t1552 = -1.*t1514;
  t1553 = 1. + t1552;
  t1547 = -1.*t1519;
  t1548 = 1. + t1547;
  t1554 = -0.0265*t1553;
  t1557 = -0.025226*t1514;
  t1558 = -0.07700600000000002*t1524;
  t1559 = t1554 + t1557 + t1558;
  t1562 = -0.2375*t1553;
  t1564 = -0.314506*t1514;
  t1566 = -0.0012740000000000008*t1524;
  t1568 = t1562 + t1564 + t1566;
  t1536 = t1519*t1514;
  t1537 = -1.*t1517*t1524;
  t1541 = t1536 + t1537;
  t1535 = Sin(var1[2]);
  t1573 = t1514*t1517;
  t1574 = t1519*t1524;
  t1575 = t1573 + t1574;
  t1508 = Cos(var1[2]);
  t1577 = -0.0695*t1548;
  t1578 = -0.0265*t1517;
  t1579 = -1.*t1517*t1559;
  t1581 = t1519*t1568;
  t1583 = t1577 + t1578 + t1579 + t1581;
  t1518 = -1.*t1514*t1517;
  t1529 = -1.*t1519*t1524;
  t1530 = t1518 + t1529;
  t1549 = -0.0265*t1548;
  t1551 = 0.0695*t1517;
  t1560 = t1519*t1559;
  t1569 = t1517*t1568;
  t1570 = t1549 + t1551 + t1560 + t1569;
  t1576 = t1570*t1575;
  t1584 = t1583*t1541;
  t1585 = t1576 + t1584;
  t1543 = -1.*t1535*t1541;
  t1598 = -1.*t1583*t1530;
  t1599 = -1.*t1570*t1541;
  t1601 = t1598 + t1599;
  t1640 = -0.0265*t1519;
  t1641 = -0.0695*t1517;
  t1642 = -1.*t1519*t1559;
  t1645 = -1.*t1517*t1568;
  t1646 = t1640 + t1641 + t1642 + t1645;
  t1631 = 0.0695*t1519;
  t1632 = t1631 + t1578 + t1579 + t1581;
  t1630 = t1583*t1530;
  t1639 = t1632*t1575;
  t1647 = t1646*t1541;
  t1648 = t1570*t1541;
  t1649 = t1630 + t1639 + t1647 + t1648;
  t1656 = -1.*t1646*t1530;
  t1657 = -1.*t1570*t1530;
  t1658 = -1.*t1632*t1541;
  t1660 = -1.*t1519*t1514;
  t1661 = t1517*t1524;
  t1662 = t1660 + t1661;
  t1666 = -1.*t1583*t1662;
  t1667 = t1656 + t1657 + t1658 + t1666;
  t1531 = t1508*t1530;
  t1544 = t1531 + t1543;
  t1610 = -1.*t1535*t1530;
  t1591 = -1.*t1535*t1575;
  t1596 = t1508*t1541;
  t1597 = t1591 + t1596;
  t1697 = t1535*t1530;
  t1698 = t1697 + t1596;
  t1719 = -0.07700600000000002*t1514;
  t1720 = t1719 + t1566;
  t1722 = -0.0012740000000000008*t1514;
  t1723 = 0.07700600000000002*t1524;
  t1724 = t1722 + t1723;
  t1728 = -1.*t1517*t1720;
  t1729 = t1519*t1724;
  t1730 = t1728 + t1729;
  t1721 = t1519*t1720;
  t1725 = t1517*t1724;
  t1726 = t1721 + t1725;
  t1671 = 0.0695*t1514;
  t1673 = t1514*t1568;
  t1675 = 0.0265*t1524;
  t1676 = t1559*t1524;
  t1677 = t1671 + t1673 + t1675 + t1676;
  t1727 = t1726*t1575;
  t1731 = t1730*t1541;
  t1732 = t1630 + t1727 + t1648 + t1731;
  t1679 = -0.0265*t1514;
  t1680 = -1.*t1514*t1559;
  t1681 = 0.0695*t1524;
  t1683 = t1568*t1524;
  t1684 = t1679 + t1680 + t1681 + t1683;
  t1734 = -1.*t1730*t1530;
  t1735 = -1.*t1726*t1541;
  t1736 = t1657 + t1734 + t1735 + t1666;
  t1689 = 0.19964*t1544*t1601;
  t1691 = t1508*t1662;
  t1692 = t1610 + t1691;
  t1693 = 0.19964*t1585*t1692;
  t1699 = 0.19964*t1698*t1601;
  t1701 = t1535*t1662;
  t1703 = t1531 + t1701;
  t1705 = 0.19964*t1585*t1703;
  t1707 = t1508*t1575;
  t1708 = t1535*t1541;
  t1709 = t1707 + t1708;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.19964*t1544*t1585 + 0.19964*t1597*t1601)*var2[0] + 0.5*(0.19964*(t1543 - 1.*t1508*t1575)*t1601 + 0.19964*t1585*(-1.*t1508*t1541 + t1610))*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[2]*(0.5*(0.19964*t1649*t1698 + t1699 + t1705 + 0.19964*t1667*t1709)*var2[0] + 0.5*(0.19964*t1544*t1649 + 0.19964*t1597*t1667 + t1689 + t1693)*var2[1] + 0.5*(0.39928*t1585*t1649 + 0.39928*t1601*t1667)*var2[2] + 0.5*(0.19964*t1649*t1677 + 0.19964*t1667*t1684)*var2[5] + 0.5*(-0.015373477840000005*t1649 - 0.0002543413600000002*t1667)*var2[6]);
  p_output1[6]=var2[2]*(0.5*(t1699 + t1705 + 0.19964*t1698*t1732 + 0.19964*t1709*t1736)*var2[0] + 0.5*(t1689 + t1693 + 0.19964*t1544*t1732 + 0.19964*t1597*t1736)*var2[1] + 0.5*(0.39928*t1585*t1732 + 0.39928*t1601*t1736)*var2[2] + 0.5*(0.19964*t1585*(0.0265*t1514 - 0.0695*t1524 + t1514*t1559 - 1.*t1524*t1568 + t1524*t1720 + t1514*t1724) + 0.19964*t1601*(t1671 + t1673 + t1675 + t1676 - 1.*t1514*t1720 + t1524*t1724) + 0.19964*t1677*t1732 + 0.19964*t1684*t1736)*var2[5] + 0.5*(-0.015373477840000005*t1732 - 0.0002543413600000002*t1736)*var2[6]);
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

#include "Ce3_vec_L5_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L5_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
