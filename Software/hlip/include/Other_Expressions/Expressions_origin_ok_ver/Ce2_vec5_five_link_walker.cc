/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:35 GMT-05:00
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
  double t1268;
  double t1252;
  double t1260;
  double t1269;
  double t1330;
  double t1267;
  double t1270;
  double t1292;
  double t1210;
  double t1341;
  double t1342;
  double t1343;
  double t1350;
  double t1351;
  double t1358;
  double t1359;
  double t1360;
  double t1361;
  double t1362;
  double t1363;
  double t1369;
  double t1329;
  double t1331;
  double t1332;
  double t1337;
  double t1338;
  double t1339;
  double t1373;
  double t1374;
  double t1375;
  double t1376;
  double t1377;
  double t1378;
  double t1393;
  double t1394;
  double t1403;
  double t1404;
  double t1405;
  double t1407;
  double t1408;
  double t1409;
  double t1413;
  double t1414;
  double t1415;
  double t1419;
  double t1420;
  double t1396;
  double t1397;
  double t1398;
  double t1370;
  double t1371;
  double t1372;
  double t1390;
  double t1391;
  double t1352;
  double t1353;
  double t1357;
  double t1365;
  double t1366;
  double t1367;
  double t1380;
  double t1381;
  double t1382;
  double t1392;
  double t1395;
  double t1399;
  double t1400;
  double t1401;
  double t1406;
  double t1410;
  double t1411;
  double t1416;
  double t1417;
  double t1418;
  double t1421;
  double t1422;
  double t1424;
  double t1425;
  double t1426;
  double t1428;
  double t1429;
  double t1430;
  double t1431;
  double t1432;
  double t1450;
  double t1451;
  double t1452;
  double t1453;
  double t1454;
  double t1412;
  double t1423;
  double t1427;
  double t1433;
  double t1434;
  double t1439;
  double t1440;
  double t1441;
  double t1442;
  double t1443;
  double t1364;
  double t1368;
  double t1379;
  double t1383;
  double t1384;
  double t1459;
  double t1460;
  double t1461;
  double t1462;
  double t1463;
  t1268 = Cos(var1[3]);
  t1252 = Cos(var1[4]);
  t1260 = Sin(var1[3]);
  t1269 = Sin(var1[4]);
  t1330 = Sin(var1[2]);
  t1267 = -1.*t1252*t1260;
  t1270 = -1.*t1268*t1269;
  t1292 = t1267 + t1270;
  t1210 = Cos(var1[2]);
  t1341 = -1.*t1252;
  t1342 = 1. + t1341;
  t1343 = 0.5*t1342;
  t1350 = 0.671885*t1252;
  t1351 = t1343 + t1350;
  t1358 = t1330*t1292;
  t1359 = t1268*t1252;
  t1360 = -1.*t1260*t1269;
  t1361 = t1359 + t1360;
  t1362 = t1210*t1361;
  t1363 = t1358 + t1362;
  t1369 = t1351*t1252;
  t1329 = t1210*t1292;
  t1331 = -1.*t1268*t1252;
  t1332 = t1260*t1269;
  t1337 = t1331 + t1332;
  t1338 = t1330*t1337;
  t1339 = t1329 + t1338;
  t1373 = t1252*t1260;
  t1374 = t1268*t1269;
  t1375 = t1373 + t1374;
  t1376 = t1210*t1375;
  t1377 = t1330*t1361;
  t1378 = t1376 + t1377;
  t1393 = -1.*t1330*t1361;
  t1394 = t1329 + t1393;
  t1403 = t1351*t1260;
  t1404 = 0.171885*t1268*t1269;
  t1405 = t1403 + t1404;
  t1407 = t1268*t1351;
  t1408 = -0.171885*t1260*t1269;
  t1409 = t1407 + t1408;
  t1413 = -0.171885*t1252*t1260;
  t1414 = -0.171885*t1268*t1269;
  t1415 = t1413 + t1414;
  t1419 = 0.171885*t1268*t1252;
  t1420 = t1419 + t1408;
  t1396 = -1.*t1330*t1292;
  t1397 = t1210*t1337;
  t1398 = t1396 + t1397;
  t1370 = Power(t1252,2);
  t1371 = -0.171885*t1370;
  t1372 = t1369 + t1371;
  t1390 = -1.*t1330*t1375;
  t1391 = t1390 + t1362;
  t1352 = t1351*t1269;
  t1353 = -0.171885*t1252*t1269;
  t1357 = t1352 + t1353;
  t1365 = -1.*t1351*t1269;
  t1366 = 0.171885*t1252*t1269;
  t1367 = t1365 + t1366;
  t1380 = Power(t1269,2);
  t1381 = 0.171885*t1380;
  t1382 = t1369 + t1381;
  t1392 = 0.85216*t1363*t1391;
  t1395 = 0.85216*t1394*t1378;
  t1399 = 0.85216*t1363*t1398;
  t1400 = 0.85216*t1394*t1339;
  t1401 = t1392 + t1395 + t1399 + t1400;
  t1406 = -1.*t1405*t1361;
  t1410 = -1.*t1292*t1409;
  t1411 = t1406 + t1410;
  t1416 = t1415*t1361;
  t1417 = t1405*t1361;
  t1418 = t1292*t1409;
  t1421 = t1375*t1420;
  t1422 = t1416 + t1417 + t1418 + t1421;
  t1424 = t1405*t1375;
  t1425 = t1361*t1409;
  t1426 = t1424 + t1425;
  t1428 = -1.*t1292*t1415;
  t1429 = -1.*t1292*t1405;
  t1430 = -1.*t1361*t1420;
  t1431 = -1.*t1409*t1337;
  t1432 = t1428 + t1429 + t1430 + t1431;
  t1450 = 0.85216*t1394*t1411;
  t1451 = 0.85216*t1394*t1422;
  t1452 = 0.85216*t1426*t1398;
  t1453 = 0.85216*t1391*t1432;
  t1454 = t1450 + t1451 + t1452 + t1453;
  t1412 = 0.85216*t1363*t1411;
  t1423 = 0.85216*t1363*t1422;
  t1427 = 0.85216*t1426*t1339;
  t1433 = 0.85216*t1378*t1432;
  t1434 = t1412 + t1423 + t1427 + t1433;
  t1439 = 0.85216*t1372*t1391;
  t1440 = 0.85216*t1357*t1394;
  t1441 = 0.85216*t1367*t1394;
  t1442 = 0.85216*t1382*t1398;
  t1443 = t1439 + t1440 + t1441 + t1442;
  t1364 = 0.85216*t1357*t1363;
  t1368 = 0.85216*t1367*t1363;
  t1379 = 0.85216*t1372*t1378;
  t1383 = 0.85216*t1382*t1339;
  t1384 = t1364 + t1368 + t1379 + t1383;
  t1459 = 0.85216*t1372*t1411;
  t1460 = 0.85216*t1367*t1426;
  t1461 = 0.85216*t1382*t1422;
  t1462 = 0.85216*t1357*t1432;
  t1463 = t1459 + t1460 + t1461 + t1462;
  p_output1[0]=var2[4]*(-0.5*(1.70432*t1339*t1363 + 1.70432*t1363*t1378)*var2[0] - 0.5*t1401*var2[1] - 0.5*t1434*var2[2] - 0.5*t1384*var2[3] - 0.0732367608*t1339*var2[4]);
  p_output1[1]=var2[4]*(-0.5*t1401*var2[0] - 0.5*(1.70432*t1391*t1394 + 1.70432*t1394*t1398)*var2[1] - 0.5*t1454*var2[2] - 0.5*t1443*var2[3] - 0.0732367608*t1398*var2[4]);
  p_output1[2]=var2[4]*(-0.5*t1434*var2[0] - 0.5*t1454*var2[1] - 0.5*(1.70432*t1422*t1426 + 1.70432*t1411*t1432)*var2[2] - 0.5*t1463*var2[3] - 0.0732367608*t1422*var2[4]);
  p_output1[3]=var2[4]*(-0.5*t1384*var2[0] - 0.5*t1443*var2[1] - 0.5*t1463*var2[2] - 0.5*(1.70432*t1357*t1372 + 1.70432*t1367*t1382)*var2[3] - 0.0732367608*t1367*var2[4]);
  p_output1[4]=(-0.0732367608*t1339*var2[0] - 0.0732367608*t1398*var2[1] - 0.0732367608*t1422*var2[2] - 0.0732367608*t1367*var2[3])*var2[4];
  p_output1[5]=0;
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

#include "Ce2_vec5_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec5_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
