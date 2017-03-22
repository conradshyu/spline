/*
 * spline.h
 *
 * Cubic Spline interpolating polynomials for free energy estimates
 * Copyright (C) 2008   Conrad Shyu (conradshyu at hotmail.com)
 * Department of Physics, University of Idaho
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Author's comments
 * -----------------
 * written by Conrad Shyu (conradshyu at hotmail.com)
 *
 * first created on December 30, 2007
 * revised on September 3, 2008
 * revised on March 6, 2014
 * revised on March 11, 2014
*/

#ifndef _SPLINE_H
#define _SPLINE_H

#include <list>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

const double TOLERANCE_LEVEL = 1e-30;

typedef struct
{
    double x;   // positions on the x-axis
    double y;   // values on the y-axis, y=f(x)
} stSPLINE;

typedef struct
{
    double x3; double x2; double x1; double x0;
} stCUBIC;

class Spline
{
public:
    Spline();
    Spline( const std::list<stSPLINE>& );
    Spline( const std::vector<double>&, const std::vector<double>& );
    ~Spline() {};

    void PrintMatrix( std::ostream& = std::cout ) const;
    bool GetEstimate( const std::string&, unsigned int ) const;

    const std::vector<double>& GetPolynomial( bool = false ) const;
    const std::list<stSPLINE>& LoadData( const std::list<stSPLINE>& );
    const std::list<stSPLINE>& LoadData( const std::vector<double>&, const std::vector<double>& );

    double DoIntegral( bool = false ) const;
    double DoQuadrature( bool = false ) const;

private:
    std::list<stSPLINE> sample;
    std::list<stCUBIC> spline;
    std::vector<double> zi;
    std::vector<double> hi;
    std::vector<double> matrix;
    std::vector<double> factor;
    unsigned int size;

    void ClearData();
    void DoGaussian();
    void DoPolynomial();        // construct the polynomial using regression
    void Exchange( double&, double& ) const;
    bool SetPivot( unsigned int );

    unsigned int SetMatrix();
    unsigned int Translate( unsigned int, unsigned int ) const;
};  // class definition for cubic spline interpolating polynomial

#endif  // _SPLINE_H
