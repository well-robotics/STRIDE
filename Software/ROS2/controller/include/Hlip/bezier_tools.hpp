/*
 * @brief Functions for computing Bezier curves and their derivatives.
 * @author Jenna Reher (jreher@caltech.edu)
 * @modifier Yuhao Huang, Yicheng Zeng
 */

#ifndef BEZIER_TOOLS_HPP
#define BEZIER_TOOLS_HPP

#include <Eigen/Dense>
#include <iostream>

namespace bezier_tools {

double singleterm_bezier(int m, int k, double s);
double bezier(const Eigen::VectorXd &coeff, double s);
void bezier(const Eigen::MatrixXd &coeffs, double s, Eigen::VectorXd &out);
double dbezier(const Eigen::VectorXd &coeff, double s);
void dbezier(const Eigen::MatrixXd &coeffs, double s, Eigen::VectorXd &out);
double d2bezier(const Eigen::VectorXd &coeff, double s);
void d2bezier(const Eigen::MatrixXd &coeffs, double s, Eigen::VectorXd &out);


}

#endif // BEZIER_TOOLS_HPP
