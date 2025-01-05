/*
 * Automatically Generated from Mathematica.
 * Fri 26 Apr 2024 13:33:47 GMT-05:00
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
  double t1316;
  double t1292;
  double t1296;
  double t1324;
  double t1298;
  double t1340;
  double t1343;
  double t1345;
  double t1346;
  double t1353;
  double t1358;
  double t1366;
  double t1369;
  double t1444;
  double t1445;
  double t1453;
  double t1454;
  double t1455;
  double t1457;
  double t1447;
  double t1448;
  double t1449;
  double t1451;
  double t1344;
  double t1354;
  double t1356;
  double t1357;
  double t1370;
  double t1399;
  double t1401;
  double t1411;
  double t1416;
  double t1417;
  double t1418;
  double t1427;
  double t1433;
  double t1434;
  double t1437;
  double t1438;
  double t1441;
  double t1442;
  double t1443;
  double t1452;
  double t1458;
  double t1459;
  double t1460;
  double t1461;
  double t1462;
  double t1463;
  double t1464;
  double t1478;
  double t1479;
  double t1473;
  double t1474;
  double t1475;
  t1316 = Cos(var1[2]);
  t1292 = Cos(var1[5]);
  t1296 = Sin(var1[2]);
  t1324 = Sin(var1[5]);
  t1298 = -1.*t1292*t1296;
  t1340 = -1.*t1316*t1324;
  t1343 = t1298 + t1340;
  t1345 = -1.*t1316*t1292;
  t1346 = t1296*t1324;
  t1353 = t1345 + t1346;
  t1358 = t1316*t1292;
  t1366 = -1.*t1296*t1324;
  t1369 = t1358 + t1366;
  t1444 = -1.*t1292;
  t1445 = 1. + t1444;
  t1453 = -0.0265*t1445;
  t1454 = -0.025413*t1292;
  t1455 = -0.08282*t1324;
  t1457 = t1453 + t1454 + t1455;
  t1447 = -0.0695*t1445;
  t1448 = -0.15232*t1292;
  t1449 = -0.0010869999999999977*t1324;
  t1451 = t1447 + t1448 + t1449;
  t1344 = -0.0007505843699999984*t1343;
  t1354 = -0.0571880382*t1353;
  t1356 = t1344 + t1354;
  t1357 = 0.5*var2[5]*t1356;
  t1370 = 1.38102*t1343*t1369;
  t1399 = 1.38102*t1343*t1353;
  t1401 = t1370 + t1399;
  t1411 = 0.5*var2[1]*t1401;
  t1416 = Power(t1343,2);
  t1417 = 0.69051*t1416;
  t1418 = t1292*t1296;
  t1427 = t1316*t1324;
  t1433 = t1418 + t1427;
  t1434 = 0.69051*t1343*t1433;
  t1437 = Power(t1369,2);
  t1438 = 0.69051*t1437;
  t1441 = 0.69051*t1369*t1353;
  t1442 = t1417 + t1434 + t1438 + t1441;
  t1443 = 0.5*var2[0]*t1442;
  t1452 = t1292*t1451;
  t1458 = t1457*t1324;
  t1459 = t1452 + t1458;
  t1460 = 0.69051*t1353*t1459;
  t1461 = -1.*t1292*t1457;
  t1462 = t1451*t1324;
  t1463 = t1461 + t1462;
  t1464 = 0.69051*t1343*t1463;
  t1478 = -0.08282*t1292;
  t1479 = t1478 + t1449;
  t1473 = -0.0010869999999999977*t1292;
  t1474 = 0.08282*t1324;
  t1475 = t1473 + t1474;
  p_output1[0]=0;
  p_output1[1]=0;
  p_output1[2]=var2[1]*(t1357 + t1411 + t1443 + 0.5*(t1460 + t1464)*var2[2]);
  p_output1[3]=0;
  p_output1[4]=0;
  p_output1[5]=var2[1]*(t1357 + t1411 + t1443 + 0.5*(t1460 + t1464 + 0.69051*t1369*(t1452 + t1458 + t1324*t1475 - 1.*t1292*t1479) + 0.69051*t1343*(-1.*t1324*t1451 + t1292*t1457 + t1292*t1475 + t1324*t1479))*var2[2]);
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

#include "Ce3_vec_L4_J2_Robot_Assembly_v3_straight_leg.hh"

namespace SymFunction
{

void Ce3_vec_L4_J2_Robot_Assembly_v3_straight_leg_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
