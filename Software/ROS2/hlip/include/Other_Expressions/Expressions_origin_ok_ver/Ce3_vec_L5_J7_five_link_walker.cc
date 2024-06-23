/*
 * Automatically Generated from Mathematica.
 * Wed 20 Sep 2023 15:41:17 GMT-05:00
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
  double t1373;
  double t1347;
  double t1361;
  double t1376;
  double t1401;
  double t1336;
  double t1413;
  double t1414;
  double t1415;
  double t1426;
  double t1437;
  double t1440;
  double t1367;
  double t1378;
  double t1384;
  double t1451;
  double t1441;
  double t1458;
  double t1466;
  double t1417;
  double t1485;
  double t1487;
  double t1488;
  double t1519;
  double t1522;
  double t1523;
  double t1524;
  double t1525;
  double t1526;
  double t1527;
  double t1530;
  double t1531;
  double t1532;
  double t1533;
  double t1534;
  double t1537;
  double t1539;
  double t1540;
  double t1542;
  double t1543;
  double t1544;
  double t1563;
  double t1564;
  double t1565;
  double t1559;
  double t1560;
  double t1484;
  double t1490;
  double t1492;
  double t1496;
  double t1498;
  double t1499;
  double t1503;
  double t1504;
  double t1505;
  double t1508;
  double t1509;
  double t1511;
  double t1513;
  double t1517;
  double t1535;
  double t1545;
  double t1549;
  double t1584;
  double t1585;
  double t1586;
  double t1580;
  double t1581;
  double t1582;
  double t1551;
  t1373 = Cos(var1[5]);
  t1347 = Cos(var1[6]);
  t1361 = Sin(var1[5]);
  t1376 = Sin(var1[6]);
  t1401 = Cos(var1[2]);
  t1336 = Sin(var1[2]);
  t1413 = t1373*t1347;
  t1414 = -1.*t1361*t1376;
  t1415 = t1413 + t1414;
  t1426 = -1.*t1347*t1361;
  t1437 = -1.*t1373*t1376;
  t1440 = t1426 + t1437;
  t1367 = t1347*t1361;
  t1378 = t1373*t1376;
  t1384 = t1367 + t1378;
  t1451 = -1.*t1336*t1415;
  t1441 = t1401*t1440;
  t1458 = t1441 + t1451;
  t1466 = -1.*t1336*t1440;
  t1417 = t1401*t1415;
  t1485 = -1.*t1373*t1347;
  t1487 = t1361*t1376;
  t1488 = t1485 + t1487;
  t1519 = -0.022659*t1347;
  t1522 = -0.007367999999999986*t1376;
  t1523 = t1519 + t1522;
  t1524 = -1.*t1361*t1523;
  t1525 = -1.*t1347;
  t1526 = 1. + t1525;
  t1527 = -0.16*t1526;
  t1530 = -0.167368*t1347;
  t1531 = 0.022659*t1376;
  t1532 = t1527 + t1530 + t1531;
  t1533 = t1373*t1532;
  t1534 = t1524 + t1533;
  t1537 = -1.*t1373*t1523;
  t1539 = -1.*t1361*t1532;
  t1540 = t1537 + t1539;
  t1542 = t1373*t1523;
  t1543 = t1361*t1532;
  t1544 = t1542 + t1543;
  t1563 = 0.022659*t1347;
  t1564 = 0.007367999999999986*t1376;
  t1565 = t1563 + t1564;
  t1559 = -0.007367999999999986*t1347;
  t1560 = t1559 + t1531;
  t1484 = 0.0033974904599999994*t1458;
  t1490 = t1401*t1488;
  t1492 = t1466 + t1490;
  t1496 = -0.0011047579199999977*t1492;
  t1498 = t1484 + t1496;
  t1499 = 0.5*var2[1]*t1498;
  t1503 = t1336*t1440;
  t1504 = t1503 + t1417;
  t1505 = 0.0033974904599999994*t1504;
  t1508 = t1336*t1488;
  t1509 = t1441 + t1508;
  t1511 = -0.0011047579199999977*t1509;
  t1513 = t1505 + t1511;
  t1517 = 0.5*var2[0]*t1513;
  t1535 = t1534*t1440;
  t1545 = t1544*t1415;
  t1549 = -1.*t1544*t1440;
  t1584 = t1373*t1565;
  t1585 = -1.*t1361*t1560;
  t1586 = t1584 + t1585;
  t1580 = t1361*t1565;
  t1581 = t1373*t1560;
  t1582 = t1580 + t1581;
  t1551 = -1.*t1534*t1488;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(0.0033974904599999994*(-1.*t1336*t1384 + t1417) - 0.0011047579199999977*t1458)*var2[0] + 0.5*(0.0033974904599999994*(-1.*t1384*t1401 + t1451) - 0.0011047579199999977*(-1.*t1401*t1415 + t1466))*var2[1])*var2[6];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(t1499 + t1517 + 0.5*(-0.0011047579199999977*(t1384*t1534 + t1535 + t1415*t1540 + t1545) + 0.0033974904599999994*(-1.*t1415*t1534 - 1.*t1440*t1540 + t1549 + t1551))*var2[2])*var2[6];
  p_output1[6]=(t1499 + t1517 + 0.5*(-0.0011047579199999977*(t1535 + t1545 + t1384*t1582 + t1415*t1586) + 0.0033974904599999994*(t1549 + t1551 - 1.*t1415*t1582 - 1.*t1440*t1586))*var2[2] + 0.5*(-0.0011047579199999977*(t1347*t1523 - 1.*t1376*t1532 + t1376*t1560 + t1347*t1565) + 0.0033974904599999994*(t1376*t1523 + t1347*t1532 - 1.*t1347*t1560 + t1376*t1565))*var2[5])*var2[6];
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

#include "Ce3_vec_L5_J7_five_link_walker.hh"

namespace SymFunction
{

void Ce3_vec_L5_J7_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
