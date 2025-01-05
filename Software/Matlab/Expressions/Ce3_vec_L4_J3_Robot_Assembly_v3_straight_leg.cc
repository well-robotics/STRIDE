/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:49 GMT-05:00
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
  double t1344;
  double t1354;
  double t1370;
  double t1399;
  double t1353;
  double t1343;
  double t1434;
  double t1437;
  double t1438;
  double t1441;
  double t1416;
  double t1417;
  double t1418;
  double t1427;
  double t1453;
  double t1454;
  double t1455;
  double t1433;
  double t1447;
  double t1448;
  double t1459;
  double t1462;
  double t1463;
  double t1491;
  double t1492;
  double t1486;
  double t1487;
  double t1488;
  double t1471;
  double t1472;
  double t1473;
  double t1474;
  double t1476;
  double t1483;
  double t1489;
  double t1490;
  double t1493;
  double t1494;
  double t1496;
  double t1497;
  double t1498;
  double t1346;
  double t1358;
  double t1369;
  double t1449;
  double t1465;
  t1344 = Cos(var1[5]);
  t1354 = Sin(var1[5]);
  t1370 = -1.*t1344;
  t1399 = 1. + t1370;
  t1353 = Sin(var1[2]);
  t1343 = Cos(var1[2]);
  t1434 = -0.0265*t1399;
  t1437 = -0.025413*t1344;
  t1438 = -0.08282*t1354;
  t1441 = t1434 + t1437 + t1438;
  t1416 = -0.0695*t1399;
  t1417 = -0.15232*t1344;
  t1418 = -0.0010869999999999977*t1354;
  t1427 = t1416 + t1417 + t1418;
  t1453 = -1.*t1344*t1353;
  t1454 = -1.*t1343*t1354;
  t1455 = t1453 + t1454;
  t1433 = t1344*t1427;
  t1447 = t1441*t1354;
  t1448 = t1433 + t1447;
  t1459 = -1.*t1344*t1441;
  t1462 = t1427*t1354;
  t1463 = t1459 + t1462;
  t1491 = -0.08282*t1344;
  t1492 = t1491 + t1418;
  t1486 = -0.0010869999999999977*t1344;
  t1487 = 0.08282*t1354;
  t1488 = t1486 + t1487;
  t1471 = 0.69051*t1455*t1448;
  t1472 = t1343*t1344;
  t1473 = -1.*t1353*t1354;
  t1474 = t1472 + t1473;
  t1476 = 0.69051*t1474*t1463;
  t1483 = t1344*t1441;
  t1489 = t1344*t1488;
  t1490 = -1.*t1427*t1354;
  t1493 = t1492*t1354;
  t1494 = t1483 + t1489 + t1490 + t1493;
  t1496 = -1.*t1344*t1492;
  t1497 = t1488*t1354;
  t1498 = t1433 + t1496 + t1447 + t1497;
  t1346 = -1.*t1343*t1344;
  t1358 = t1353*t1354;
  t1369 = t1346 + t1358;
  t1449 = 0.69051*t1369*t1448;
  t1465 = 0.69051*t1455*t1463;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(0.5*(t1471 + t1476)*var2[0] + 0.5*(t1449 + t1465)*var2[1])*var2[2];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[2]*(0.5*(t1471 + t1476 + 0.69051*t1474*t1494 + 0.69051*(t1344*t1353 + t1343*t1354)*t1498)*var2[0] + 0.5*(t1449 + t1465 + 0.69051*t1455*t1494 + 0.69051*t1474*t1498)*var2[1] + 0.5*(1.38102*t1448*t1494 + 1.38102*t1463*t1498)*var2[2] + 0.5*(-0.0571880382*t1494 - 0.0007505843699999984*t1498)*var2[5]);
  p_output1[6]=0;
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

#include "Ce3_vec_L4_J3_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L4_J3_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
