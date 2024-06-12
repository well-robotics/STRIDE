/*
 * MIT License
 * 
 * Copyright (c) 2020 Jenna Reher (jreher@caltech.edu)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
*/
#include <bezier_tools.hpp>

namespace bezier_tools {

double factorial(int n)
{
    double out = 1.;
    for (int i = 2; i <= n; ++i)
        out *= i;
    return out;
}

double nchoosek(int n, int k)
{
    return factorial(n) / (factorial(k) * factorial(n - k));
}

double singleterm_bezier(int m, int k, double s)
{
    if (k == 0)
        return nchoosek(m, k) * std::pow(1 - s, m - k);
    else if (k == m)
        return nchoosek(m, k) * std::pow(s, k);
    else
        return nchoosek(m, k) * std::pow(s, k) * std::pow(1 - s, m - k);
}

double bezier(const Eigen::VectorXd &coeff, double s)
{
    int m = coeff.size() - 1;
    double fcn = 0.;
    for (int k = 0; k <= m; ++k)
    {
        fcn += coeff(k) * singleterm_bezier(m, k, s);
    }
    return fcn;
}

void bezier(const Eigen::MatrixXd &coeffs, double s, Eigen::VectorXd &out)
{
    for (int i = 0; i < coeffs.rows(); ++i)
        out(i) = bezier(coeffs.row(i), s);
}

double diff_coeff(const Eigen::VectorXd &coeff, Eigen::VectorXd &dcoeff)
{
    int m = coeff.size() - 1;
    Eigen::MatrixXd A(m, m + 1);
    A.setZero();

    for (int i = 0; i < m; ++i)
    {
        A(i, i) = -(m - i) * nchoosek(m, i) / nchoosek(m - 1, i);
        A(i, i + 1) = (i + 1) * nchoosek(m, i + 1) / nchoosek(m - 1, i);
    }
    A(m - 1, m) = m * nchoosek(m, m);
    dcoeff = A * coeff;
    return dcoeff(0);//added March 4th
}

double dbezier(const Eigen::VectorXd &coeff, double s)
{
    Eigen::VectorXd dcoeff;
    dcoeff.resizeLike(coeff);
    // std::cout << "dcoeff size" << dcoeff.size() << std::endl;
    diff_coeff(coeff, dcoeff);
    return bezier(dcoeff, s);
}

void dbezier(const Eigen::MatrixXd &coeffs, double s, Eigen::VectorXd &out)
{
    for (int i = 0; i < coeffs.rows(); ++i)
        out(i) = dbezier(coeffs.row(i), s);
}

double d2bezier(const Eigen::VectorXd &coeff, double s)
{
    Eigen::VectorXd dcoeff, d2coeff;
    dcoeff.resizeLike(coeff);
    d2coeff.resizeLike(coeff);
    diff_coeff(coeff, dcoeff);
    diff_coeff(dcoeff, d2coeff);
    return bezier(d2coeff, s);
}

void d2bezier(const Eigen::MatrixXd &coeffs, double s, Eigen::VectorXd &out)
{
    for (int i = 0; i < coeffs.rows(); ++i)
        out(i) = d2bezier(coeffs.row(i), s);
}

}











