/*
 * Automatically Generated from Mathematica.
 * Mon 24 Jul 2023 21:28:34 GMT-05:00
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
  double t1099;
  double t1051;
  double t1065;
  double t1102;
  double t1121;
  double t1010;
  double t1091;
  double t1106;
  double t1114;
  double t1223;
  double t1224;
  double t1225;
  double t1226;
  double t1227;
  double t1117;
  double t1122;
  double t1123;
  double t1165;
  double t1208;
  double t1209;
  double t1211;
  double t1216;
  double t1217;
  double t1253;
  double t1254;
  double t1255;
  double t1231;
  double t1240;
  double t1241;
  double t1242;
  double t1243;
  double t1244;
  double t1257;
  double t1258;
  double t1259;
  double t1261;
  double t1262;
  double t1263;
  double t1264;
  double t1265;
  double t1266;
  double t1283;
  double t1284;
  double t1299;
  double t1300;
  double t1301;
  double t1303;
  double t1304;
  double t1305;
  double t1309;
  double t1310;
  double t1311;
  double t1286;
  double t1287;
  double t1288;
  double t1276;
  double t1277;
  double t1278;
  double t1228;
  double t1229;
  double t1230;
  double t1246;
  double t1247;
  double t1248;
  double t1249;
  double t1256;
  double t1280;
  double t1281;
  double t1271;
  double t1272;
  double t1273;
  double t1274;
  double t1275;
  double t1279;
  double t1282;
  double t1285;
  double t1289;
  double t1290;
  double t1291;
  double t1293;
  double t1294;
  double t1295;
  double t1296;
  double t1297;
  double t1302;
  double t1306;
  double t1307;
  double t1312;
  double t1313;
  double t1314;
  double t1315;
  double t1316;
  double t1318;
  double t1319;
  double t1320;
  double t1322;
  double t1323;
  double t1324;
  double t1325;
  double t1326;
  double t1344;
  double t1345;
  double t1346;
  double t1347;
  double t1348;
  double t1349;
  double t1298;
  double t1308;
  double t1317;
  double t1321;
  double t1327;
  double t1328;
  double t1333;
  double t1334;
  double t1335;
  double t1336;
  double t1222;
  double t1245;
  double t1250;
  double t1251;
  double t1354;
  double t1355;
  double t1356;
  t1099 = Cos(var1[3]);
  t1051 = Cos(var1[4]);
  t1065 = Sin(var1[3]);
  t1102 = Sin(var1[4]);
  t1121 = Sin(var1[2]);
  t1010 = Cos(var1[2]);
  t1091 = -1.*t1051*t1065;
  t1106 = -1.*t1099*t1102;
  t1114 = t1091 + t1106;
  t1223 = -1.*t1051;
  t1224 = 1. + t1223;
  t1225 = 0.5*t1224;
  t1226 = 0.671885*t1051;
  t1227 = t1225 + t1226;
  t1117 = t1010*t1114;
  t1122 = -1.*t1099*t1051;
  t1123 = t1065*t1102;
  t1165 = t1122 + t1123;
  t1208 = t1121*t1165;
  t1209 = t1117 + t1208;
  t1211 = -1.*t1099*t1121;
  t1216 = -1.*t1010*t1065;
  t1217 = t1211 + t1216;
  t1253 = t1010*t1099;
  t1254 = -1.*t1121*t1065;
  t1255 = t1253 + t1254;
  t1231 = t1121*t1114;
  t1240 = t1099*t1051;
  t1241 = -1.*t1065*t1102;
  t1242 = t1240 + t1241;
  t1243 = t1010*t1242;
  t1244 = t1231 + t1243;
  t1257 = t1099*t1121;
  t1258 = t1010*t1065;
  t1259 = t1257 + t1258;
  t1261 = t1051*t1065;
  t1262 = t1099*t1102;
  t1263 = t1261 + t1262;
  t1264 = t1010*t1263;
  t1265 = t1121*t1242;
  t1266 = t1264 + t1265;
  t1283 = -1.*t1121*t1242;
  t1284 = t1117 + t1283;
  t1299 = t1227*t1065;
  t1300 = 0.171885*t1099*t1102;
  t1301 = t1299 + t1300;
  t1303 = t1099*t1227;
  t1304 = -0.171885*t1065*t1102;
  t1305 = t1303 + t1304;
  t1309 = -1.*t1227*t1065;
  t1310 = -0.171885*t1099*t1102;
  t1311 = t1309 + t1310;
  t1286 = -1.*t1121*t1114;
  t1287 = t1010*t1165;
  t1288 = t1286 + t1287;
  t1276 = -1.*t1010*t1099;
  t1277 = t1121*t1065;
  t1278 = t1276 + t1277;
  t1228 = t1227*t1102;
  t1229 = -0.171885*t1051*t1102;
  t1230 = t1228 + t1229;
  t1246 = t1227*t1051;
  t1247 = Power(t1102,2);
  t1248 = 0.171885*t1247;
  t1249 = t1246 + t1248;
  t1256 = 6.8522*t1217*t1255;
  t1280 = -1.*t1121*t1263;
  t1281 = t1280 + t1243;
  t1271 = Power(t1217,2);
  t1272 = 3.4261*t1271;
  t1273 = 3.4261*t1217*t1259;
  t1274 = Power(t1255,2);
  t1275 = 3.4261*t1274;
  t1279 = 3.4261*t1255*t1278;
  t1282 = 0.85216*t1244*t1281;
  t1285 = 0.85216*t1284*t1266;
  t1289 = 0.85216*t1244*t1288;
  t1290 = 0.85216*t1284*t1209;
  t1291 = t1272 + t1273 + t1275 + t1279 + t1282 + t1285 + t1289 + t1290;
  t1293 = Power(t1099,2);
  t1294 = 0.1494*t1293;
  t1295 = Power(t1065,2);
  t1296 = 0.1494*t1295;
  t1297 = t1294 + t1296;
  t1302 = -1.*t1301*t1242;
  t1306 = -1.*t1114*t1305;
  t1307 = t1302 + t1306;
  t1312 = t1311*t1242;
  t1313 = t1301*t1242;
  t1314 = t1114*t1305;
  t1315 = t1263*t1305;
  t1316 = t1312 + t1313 + t1314 + t1315;
  t1318 = t1301*t1263;
  t1319 = t1242*t1305;
  t1320 = t1318 + t1319;
  t1322 = -1.*t1114*t1311;
  t1323 = -1.*t1114*t1301;
  t1324 = -1.*t1242*t1305;
  t1325 = -1.*t1305*t1165;
  t1326 = t1322 + t1323 + t1324 + t1325;
  t1344 = 3.4261*t1278*t1297;
  t1345 = 0.85216*t1284*t1307;
  t1346 = 0.85216*t1284*t1316;
  t1347 = 0.85216*t1320*t1288;
  t1348 = 0.85216*t1281*t1326;
  t1349 = t1344 + t1345 + t1346 + t1347 + t1348;
  t1298 = 3.4261*t1217*t1297;
  t1308 = 0.85216*t1244*t1307;
  t1317 = 0.85216*t1244*t1316;
  t1321 = 0.85216*t1320*t1209;
  t1327 = 0.85216*t1266*t1326;
  t1328 = t1298 + t1308 + t1317 + t1321 + t1327;
  t1333 = 0.51185934*t1278;
  t1334 = 0.85216*t1230*t1284;
  t1335 = 0.85216*t1249*t1288;
  t1336 = t1333 + t1334 + t1335;
  t1222 = 0.51185934*t1217;
  t1245 = 0.85216*t1230*t1244;
  t1250 = 0.85216*t1249*t1209;
  t1251 = t1222 + t1245 + t1250;
  t1354 = 0.85216*t1249*t1316;
  t1355 = 0.85216*t1230*t1326;
  t1356 = t1354 + t1355;
  p_output1[0]=var2[3]*(-0.5*(1.70432*t1209*t1244 + t1256 + 6.8522*t1255*t1259 + 1.70432*t1244*t1266)*var2[0] - 0.5*t1291*var2[1] - 0.5*t1328*var2[2] - 0.5*t1251*var2[3] - 0.0732367608*t1209*var2[4]);
  p_output1[1]=var2[3]*(-0.5*t1291*var2[0] - 0.5*(t1256 + 6.8522*t1217*t1278 + 1.70432*t1281*t1284 + 1.70432*t1284*t1288)*var2[1] - 0.5*t1349*var2[2] - 0.5*t1336*var2[3] - 0.0732367608*t1288*var2[4]);
  p_output1[2]=var2[3]*(-0.5*t1328*var2[0] - 0.5*t1349*var2[1] - 0.5*(1.70432*t1316*t1320 + 1.70432*t1307*t1326)*var2[2] - 0.5*t1356*var2[3] - 0.0732367608*t1316*var2[4]);
  p_output1[3]=(-0.5*t1251*var2[0] - 0.5*t1336*var2[1] - 0.5*t1356*var2[2])*var2[3];
  p_output1[4]=(-0.0732367608*t1209*var2[0] - 0.0732367608*t1288*var2[1] - 0.0732367608*t1316*var2[2])*var2[3];
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

#include "Ce2_vec4_five_link_walker.hh"

namespace SymFunction
{

void Ce2_vec4_five_link_walker_raw(double *p_output1, const double *var1,const double *var2)
{
  // Call Subroutines
  output1(p_output1, var1, var2);

}

}

#endif // MATLAB_MEX_FILE
