/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:41:32 GMT-05:00
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
  double t1391;
  double t1367;
  double t1441;
  double t1459;
  double t1458;
  double t1460;
  double t1464;
  double t1469;
  double t1470;
  double t1471;
  double t1480;
  double t1487;
  double t1527;
  double t1531;
  double t1530;
  double t1533;
  double t1535;
  double t1541;
  double t1542;
  double t1543;
  double t1548;
  double t1554;
  double t1468;
  double t1472;
  double t1476;
  double t1479;
  double t1482;
  double t1483;
  double t1484;
  double t1490;
  double t1492;
  double t1496;
  double t1503;
  double t1504;
  double t1505;
  double t1508;
  double t1509;
  double t1511;
  double t1519;
  double t1522;
  double t1574;
  double t1575;
  double t1576;
  double t1536;
  double t1545;
  double t1546;
  double t1547;
  double t1550;
  double t1552;
  double t1553;
  double t1555;
  double t1556;
  double t1557;
  double t1558;
  double t1559;
  double t1561;
  double t1562;
  double t1563;
  double t1564;
  double t1566;
  double t1567;
  double t1591;
  double t1592;
  double t1593;
  t1391 = Sin(var1[2]);
  t1367 = Cos(var1[2]);
  t1441 = Cos(var1[3]);
  t1459 = Sin(var1[3]);
  t1458 = -1.*t1441*t1391;
  t1460 = -1.*t1367*t1459;
  t1464 = t1458 + t1460;
  t1469 = -1.*t1367*t1441;
  t1470 = t1391*t1459;
  t1471 = t1469 + t1470;
  t1480 = Cos(var1[4]);
  t1487 = Sin(var1[4]);
  t1527 = Cos(var1[5]);
  t1531 = Sin(var1[5]);
  t1530 = -1.*t1527*t1391;
  t1533 = -1.*t1367*t1531;
  t1535 = t1530 + t1533;
  t1541 = -1.*t1367*t1527;
  t1542 = t1391*t1531;
  t1543 = t1541 + t1542;
  t1548 = Cos(var1[6]);
  t1554 = Sin(var1[6]);
  t1468 = -0.078722*t1464;
  t1472 = -0.00102*t1471;
  t1476 = t1468 + t1472;
  t1479 = -6.3425574000000005*t1476;
  t1482 = -1.*t1480;
  t1483 = 1. + t1482;
  t1484 = -0.16*t1483*t1464;
  t1490 = 0.16*t1471*t1487;
  t1492 = t1480*t1471;
  t1496 = -1.*t1464*t1487;
  t1503 = t1492 + t1496;
  t1504 = -0.022663*t1503;
  t1505 = t1480*t1464;
  t1508 = t1471*t1487;
  t1509 = t1505 + t1508;
  t1511 = -0.167371*t1509;
  t1519 = t1484 + t1490 + t1504 + t1511;
  t1522 = -1.4709113999999999*t1519;
  t1574 = t1367*t1441;
  t1575 = -1.*t1391*t1459;
  t1576 = t1574 + t1575;
  t1536 = -0.078865*t1535;
  t1545 = -0.001112*t1543;
  t1546 = t1536 + t1545;
  t1547 = -6.3204849*t1546;
  t1550 = -1.*t1548;
  t1552 = 1. + t1550;
  t1553 = -0.16*t1552*t1535;
  t1555 = 0.16*t1543*t1554;
  t1556 = t1548*t1543;
  t1557 = -1.*t1535*t1554;
  t1558 = t1556 + t1557;
  t1559 = -0.022659*t1558;
  t1561 = t1548*t1535;
  t1562 = t1543*t1554;
  t1563 = t1561 + t1562;
  t1564 = -0.167368*t1563;
  t1566 = t1553 + t1555 + t1559 + t1564;
  t1567 = -1.4709113999999999*t1566;
  t1591 = t1367*t1527;
  t1592 = -1.*t1391*t1531;
  t1593 = t1591 + t1592;
  p_output1[0]=0;
  p_output1[1]=-27.864422100000002;
  p_output1[2]=-12.259557000000001*(-0.000042*t1367 - 0.007972*t1391) + t1479 + t1522 + t1547 + t1567;
  p_output1[3]=t1479 + t1522;
  p_output1[4]=-1.4709113999999999*(0.16*t1464*t1480 - 0.16*t1487*t1576 - 0.022663*(t1496 - 1.*t1480*t1576) - 0.167371*(t1505 - 1.*t1487*t1576));
  p_output1[5]=t1547 + t1567;
  p_output1[6]=-1.4709113999999999*(0.16*t1535*t1548 - 0.16*t1554*t1593 - 0.022659*(t1557 - 1.*t1548*t1593) - 0.167368*(t1561 - 1.*t1554*t1593));
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

#include "Ge_vec_five_link_walker.hh"

namespace SymFunction
{

void Ge_vec_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
