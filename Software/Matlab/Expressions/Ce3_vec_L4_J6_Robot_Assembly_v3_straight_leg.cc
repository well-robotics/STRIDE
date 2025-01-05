/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:56 GMT-05:00
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
  double t1417;
  double t1358;
  double t1369;
  double t1418;
  double t1416;
  double t1433;
  double t1434;
  double t1437;
  double t1438;
  double t1447;
  double t1448;
  double t1455;
  double t1462;
  double t1463;
  double t1467;
  double t1468;
  double t1472;
  double t1474;
  double t1477;
  double t1478;
  double t1480;
  double t1486;
  double t1487;
  double t1502;
  double t1500;
  double t1501;
  double t1503;
  double t1505;
  double t1506;
  double t1488;
  double t1489;
  double t1491;
  double t1493;
  double t1495;
  double t1497;
  double t1498;
  t1417 = Cos(var1[2]);
  t1358 = Cos(var1[5]);
  t1369 = Sin(var1[2]);
  t1418 = Sin(var1[5]);
  t1416 = -1.*t1358*t1369;
  t1433 = -1.*t1417*t1418;
  t1434 = t1416 + t1433;
  t1437 = -0.0571880382*t1434;
  t1438 = t1417*t1358;
  t1447 = -1.*t1369*t1418;
  t1448 = t1438 + t1447;
  t1455 = -0.0007505843699999984*t1448;
  t1462 = t1437 + t1455;
  t1463 = 0.5*var2[0]*t1462;
  t1467 = -0.0007505843699999984*t1434;
  t1468 = -1.*t1417*t1358;
  t1472 = t1369*t1418;
  t1474 = t1468 + t1472;
  t1477 = -0.0571880382*t1474;
  t1478 = t1467 + t1477;
  t1480 = 0.5*var2[1]*t1478;
  t1486 = -1.*t1358;
  t1487 = 1. + t1486;
  t1502 = -0.0010869999999999977*t1418;
  t1500 = -0.0695*t1487;
  t1501 = -0.15232*t1358;
  t1503 = t1500 + t1501 + t1502;
  t1505 = -0.08282*t1358;
  t1506 = t1505 + t1502;
  t1488 = -0.0265*t1487;
  t1489 = -0.025413*t1358;
  t1491 = -0.08282*t1418;
  t1493 = t1488 + t1489 + t1491;
  t1495 = -0.0010869999999999977*t1358;
  t1497 = 0.08282*t1418;
  t1498 = t1495 + t1497;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=(t1463 + t1480)*var2[5];
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=(t1463 + t1480 + 0.5*(-0.0007505843699999984*(t1418*t1493 + t1418*t1498 + t1358*t1503 - 1.*t1358*t1506) - 0.0571880382*(t1358*t1493 + t1358*t1498 - 1.*t1418*t1503 + t1418*t1506))*var2[2])*var2[5];
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

#include "Ce3_vec_L4_J6_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L4_J6_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
