/*
 * spline.cpp
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
*/

#include <spline.h>

/*
 * default class constructor
*/
Spline::Spline()
{
    ClearData();
}   // end of class constructor

/*
 * class constructor
*/
Spline::Spline(
    const std::list<stSPLINE>& _sample )
{
    LoadData( _sample );
}   // end of class constructor

/*
 * class constructor
*/
Spline::Spline(
    const std::vector<double>& _x,
    const std::vector<double>& _y )
{
    LoadData( _x, _y );
}   // end of class constructor

double Spline::DoIntegral(
    bool _print ) const
{
    double area = 0.0;
    std::list<stSPLINE>::const_iterator a = sample.begin();
    std::list<stSPLINE>::const_iterator b = sample.begin(); b++;

    for ( std::list<stCUBIC>::const_iterator i = spline.begin(); !( i == spline.end() ); i++ )
    {
        area += ( *i ).x3 * ( pow( ( *b ).x, 4.0 ) - pow( ( *a ).x, 4.0 ) ) * 0.25 +
            ( *i ).x2 * ( pow( ( *b ).x, 3.0 ) - pow( ( *a ).x, 3.0 ) ) * ( 1.0 / 3.0 ) +
            ( *i ).x1 * ( pow( ( *b ).x, 2.0 ) - pow( ( *a ).x, 2.0 ) ) * 0.5 +
            ( *i ).x0 * ( ( *b ).x - ( *a ).x );
        a++; b++;
    }   // process each polynomial

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the result

    return( area );
}   // end of DoIntegral()

/*
 * calculate the area under the curve using quadrature
*/
double Spline::DoQuadrature(
    bool _print ) const
{
    std::list<stSPLINE>::const_iterator a = sample.begin();
    std::list<stSPLINE>::const_iterator b = sample.begin(); b++;
    double area = 0.0;

    while ( !( b == sample.end() ) )
    {
        area += ( ( *b ).y + ( *a ).y ) * 0.5 * ( ( *b ).x - ( *a ).x );
        a = b; b++;
    }   // iterate through the entire list

    if ( _print )
    {
        printf( "area under the curve: %.8f\n", area );
    }   // print out the integration result

    return( area );
}   // end of DoQuadrature()

/*
 * construct the cubic spline polynomials
*/
void Spline::DoPolynomial()
{
    SetMatrix(); DoGaussian();  // construct the matrix and solve by gaussina elimination

    stCUBIC cubic; spline.clear(); factor.clear();
    std::list<stSPLINE>::iterator a = sample.begin();
    std::list<stSPLINE>::iterator b = sample.begin(); b++;

    for ( unsigned int i = 1; i < zi.size(); a++, b++, ++i )
    {
        cubic.x3 = ( 1.0 / ( 6.0 * hi[ i - 1 ] ) ) * ( zi[ i ] - zi[ i - 1 ] );
        cubic.x2 = ( 1.0 / ( 2.0 * hi[ i - 1 ] ) ) * ( zi[ i - 1 ] * ( *b ).x - zi[ i ] * ( *a ).x );
        cubic.x1 = ( 1.0 / ( 2.0 * hi[ i - 1 ] ) ) *
            ( zi[ i ] * pow( ( *a ).x, 2.0 ) - zi[ i - 1 ] * pow( ( *b ).x, 2.0 ) ) +
            ( 1.0 / hi[ i - 1 ] ) * ( ( *b ).y - ( *a ).y ) -
            ( hi[ i - 1 ] / 6.0 ) * ( zi[ i ] - zi[ i - 1 ] );
        cubic.x0 = ( 1.0 / ( 6.0 * hi[ i - 1 ] ) ) *
            ( zi[ i - 1 ] * pow( ( *b ).x, 3.0 ) - zi[ i ] * pow( ( *a ).x, 3.0 ) ) +
            ( 1.0 / hi[ i - 1 ] ) * ( ( *a ).y * ( *b ).x - ( *b ).y * ( *a ).x ) +
            ( hi[ i - 1 ] / 6.0 ) * ( zi[ i ] * ( *a ).x - zi[ i - 1 ] * ( *b ).x );

        spline.push_back( cubic );
    }   // construct the cubic spline polynomial for each interval

    for ( std::list<stCUBIC>::iterator j = spline.begin(); !( j == spline.end() ); j++ )
    {
        factor.push_back( ( *j ).x0 ); factor.push_back( ( *j ).x1 );
        factor.push_back( ( *j ).x2 ); factor.push_back( ( *j ).x3 );
    }   // save a copy of the polynomial coefficients
}   // end of DoPolynomial()

/*
 * perform the gaussian elimination
*/
void Spline::DoGaussian()
{
    double ratio; unsigned int u;

    for ( unsigned int s = 0; s < size; ++s )
    {
        if ( !SetPivot( s ) )
        {
            std::cout << "Error: matrix is singular" << std::endl; exit( 1 );
        }   // apply partial pivoting to the matrix and check for singularity

        for ( unsigned int i = ( s + 1 ); i < size; ++i )
        {
            ratio = matrix[ Translate( i, s ) ] / matrix[ Translate( s, s ) ];
            zi[ i + 1 ] -= zi[ s + 1 ] * ratio;

            for ( unsigned int j = s; j < size; ++j )
            {
                matrix[ Translate( i, j ) ] -= ( matrix[ Translate( s, j ) ] * ratio );
            }   // successively remove previous terms
        }   // perform forward elimination on the matrix
    }   // perform gauss elimination with partial pivoting

    for ( unsigned int offset = 0; offset < size; ++offset )
    {
        u = size - offset - 1;
        ratio = zi[ u + 1 ] / matrix[ Translate( u, u ) ];
        matrix[ Translate( u, u ) ] = 1.0; zi[ u + 1 ] = ratio;

        for ( unsigned int v = 0; v < u; ++v )
        {
            zi[ v + 1 ] -= ( matrix[ Translate( v, u ) ] * ratio );
            matrix[ Translate( v, u ) ] = 0.0;
        }   // update the solution to the linear equations
    }   // perform backward substitution
}   // end of DoGaussian()

/*
 * construct the matrix for calculating the z values
 *
 * | 2*(h0+h1)    h1         0         0     || z1 |     | (y2-y1)/h1 - (y1-y0)/h0 |
 * |    h1     2*(h1+h2)    h2         0     || z2 | = 6*| (y3-y2)/h2 - (y2-y1)/h1 |
 * |     0        h2     2*(h2+h3)    h3     || z3 |     | (y4-y3)/h3 - (y3-y2)/h2 |
 * |     0         0        h3     2*(h3+h4) || z4 |     | (y5-y4)/h4 - (y4-y3)/h3 |
 *
 * ------------ variable: matrix -------------           ------ variable: zi -------
 *
 * for natural spline: z0 = z5 = 0
 * f(xi)=yi, and hi is the size of the interval
 *
 * note: the construction of the linear equation matrix has been verified to produce
 * correct results on May 12, 2008
*/
unsigned int Spline::SetMatrix()
{
    size = sample.size() - 2;       // ignoare z0 and zn

    matrix.clear(); matrix.resize( size * size, 0.0 );
    zi.clear(); zi.resize( sample.size(), 0.0 );
    hi.clear(); hi.resize( sample.size() - 1, 0.0 );

    std::vector<double> yi( sample.size() - 1, 0.0 );

    std::list<stSPLINE>::iterator ha = sample.begin();
    std::list<stSPLINE>::iterator hb = sample.begin(); hb++;

    for ( unsigned int h = 0; !( hb == sample.end() ); ha++, hb++, ++h )
    {
        hi[ h ] = ( *hb ).x - ( *ha ).x;
        yi[ h ] = ( ( *hb ).y - ( *ha ).y ) / hi[ h ];
    }   // calculate the difference between two adjacent intervals

    for ( unsigned int i = 0; i < size; ++i )
    {
        matrix[ Translate( i, i ) ] = 2.0 * ( hi[ i ] + hi[ i + 1 ] );
        zi[ i + 1 ] = 6.0 * ( yi[ i + 1 ] - yi[ i ] );
    }   // assign the diagonal elements and construct the matrix for z_i

    for ( unsigned int j = 1; j < size; ++j )
    {
        matrix[ Translate( j - 1, j ) ] = hi[ j ];  // upper diagonal elements
        matrix[ Translate( j, j - 1 ) ] = hi[ j ];  // lower diagonal elements
    }   // assign the upper and lower diagonal elements

    return( matrix.size() );
}   // end of SetMatrix()

/*
 * set the matrix pivoting elements
 * note: matrix pivoting has been verified to work correctly on december 29, 2007
*/
bool Spline::SetPivot(
    unsigned int _r )       // current diagonal position
{
    for ( unsigned i = ( _r + 1 ); i < size; ++i )
    {
        if ( !( matrix[ Translate( i, _r ) ] > matrix[ Translate( _r, _r ) ] ) )
        {
            continue;
        }   // search the largest value for pivot

        for ( unsigned int j = 0; j < size; ++j )
        {
            Exchange( matrix[ Translate( i, j ) ], matrix[ Translate( _r, j ) ] );
        }   // swap the elements in the matrix

        Exchange( zi[ i ], zi[ _r ] );
    }   // perform partial pivoting on the matrix

    return( ( fabs( matrix[ Translate( _r, _r ) ] ) > TOLERANCE_LEVEL ) ? true : false );
}   // end of SetPivot()

/*
 * get the estimate from the polynomial
*/
bool Spline::GetEstimate(
    const std::string& _file,
    const unsigned int _step ) const
{
    std::ofstream ofs( _file.c_str(), std::ios::trunc );

    if ( ofs.bad() )
    {
        std::cout << "file " << _file << "cannot be opened" << std::endl;
        return( false );
    }   // make sure the file stream has been opened successfully

    std::list<stSPLINE>::const_iterator j; std::list<stCUBIC>::const_iterator i;
    char buffer[ 80 ]; double y = 0.0; double x = 0.0;
    double step = 1.0 / static_cast<double>( _step );

    for ( unsigned int s = 0; !( s > _step ); ++s )
    {
        i = spline.begin(); j = sample.begin(); j++;    // reset the record pointers

        while ( x > ( ( *j ).x + step ) )
        {
            i++; j++;
        }   // search for the correct interval

        y = ( *i ).x3 * pow( x, 3.0 ) + ( *i ).x2 * pow( x, 2.0 ) + ( *i ).x1 * x + ( *i ).x0;
        sprintf( buffer, "%.4f, %.8f", x, y );
        ofs << buffer << std::endl; y = 0.0; x += step;
    }   // iterate through the entire interval

    ofs.close(); return( true );
}   // end of GetEstimate()

/*
 * reset and initialize essential variables
*/
const std::list<stSPLINE>& Spline::LoadData(
    const std::list<stSPLINE>& _sample )
{
    stSPLINE unit; ClearData();

    for ( std::list<stSPLINE>::const_iterator i = _sample.begin(); !( i == _sample.end() ); i++ )
    {
        unit.x = ( *i ).x; unit.y = ( *i ).y; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation using cubic spline
    DoPolynomial(); return( sample );
}   // end of LoadData()

/*
 * reset and initialize essential variables
*/
const std::list<stSPLINE>& Spline::LoadData(
    const std::vector<double>& _x,
    const std::vector<double>& _y )
{
    stSPLINE unit; ClearData();

    for ( unsigned int i = 0; i < _x.size(); ++i )
    {
        unit.x = _x[ i ]; unit.y = _y[ i ]; sample.push_back( unit );
    }   // save a local copy of the data

    // perform interpolation using cubic spline
    DoPolynomial(); return( sample );
}   // end of LoadData()

/*
 * print out the cubic spline polynomials
*/
const std::vector<double>& Spline::GetPolynomial(
    bool _print ) const
{
    std::list<stCUBIC>::const_iterator i;
    std::list<stSPLINE>::const_iterator s;
    double r0, r1;

    if ( _print )
    {
        printf( "   Interval, Polynomial coefficients\n" );

        for ( i = spline.begin(), s = sample.begin(); !( i == spline.end() ); i++ )
        {
            r0 = ( *s ).x; s++; r1 = ( *s ).x;
            printf( "%1.2f - %1.2f, %.8f %.8f %.8f %.8f\n",
                r0, r1, ( *i ).x0, ( *i ).x1, ( *i ).x2, ( *i ).x3 );
        }   // format the cubic spline polynomials
    }   // print out the cubic spline polynomials

    return( factor );
}   // end of PrintPolynomial()

/*
 * print out the matrix
*/
void Spline::PrintMatrix(
    std::ostream& _os ) const
{
    char buffer[ 80 ];

    for ( unsigned int i = 0; i < size; ++i )
    {
        for ( unsigned int j = 0; j < size; ++j )
        {
            sprintf( buffer, "%.4f ", matrix[ Translate( i, j ) ] ); _os << buffer;
        }   // stupid c++ iostream, can't format the output easily

        sprintf( buffer, "| %.6f", zi[ i + 1 ] ); _os << buffer << std::endl;
    }   // print the matrix
}   // end of PrintMatrix()

/*
 * clear all contents
*/
void Spline::ClearData()
{
    sample.clear(); zi.clear(); hi.clear();
}   // end of ClearData()

/*
 * swap the contents of two variables
*/
void Spline::Exchange(
    double& _a, double& _b ) const
{
    double swap = _a; _a = _b; _b = swap;
}   // end of Exchange()

/*
 * translate two dimensional coordinate into one dimension
*/
unsigned int Spline::Translate(
    unsigned int _r, unsigned int _c ) const
{
    return( _r * size + _c );
}   // end of Translate()
