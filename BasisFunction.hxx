////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _BASISFUNCTION_HXX
#define _BASISFUNCTION_HXX 1

#include <cmath>
#include <vector>
#include <iostream>
using namespace std;

#include <Point3.h>
#include <Vector3.h>
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#define LINEAR    1
#define QUADRATIC 2
#define CUBIC     3

#define LINEAR_COEFF     4
#define QUADRATIC_COEFF 10
#define CUBIC_COEFF     20

#define LINEAR_COEFF_ALL     8
#define QUADRATIC_COEFF_ALL 14
#define CUBIC_COEFF_ALL     24

template <typename T>
class BasisFunction {

public:

  BasisFunction(){};
  virtual ~BasisFunction(){};

  // p は相対座標
  T poly( Point3<T>& p ) const {
    return poly_( p.x, p.y, p.z );
  };
  
  void polyG( Point3<T>& p, Vector3<T>& g ) {
    
    T ga[3];
    polyG_( ga, p.x, p.y, p.z );
    g.set( ga[0], ga[1], ga[2] );

  };

  void resize( unsigned int n ) { c_.resize( n ); };
  void clear() { c_.clear(); };

  void setCoeff( int i, T f ) { c_[i] = f; };
  T coeff( int i ) const { return c_[i]; };
  std::vector<T>& coeffs() { return c_; };

private:

  std::vector<T> c_;

  virtual T poly_( T vx, T vy, T vz ) const { return .0f; };
  virtual void polyG_( T g[3], T vx, T vy, T vz ) {};


};

template <typename T>
class Cubic : public BasisFunction<T> {

public:

  Cubic(){
    resize( CUBIC_COEFF );
  };

  ~Cubic(){
    clear();
  };
  
//   std::vector<T>& coeff() { return c_; };

  void setCoeffs( std::vector<T>& c ) {
    setCoeffs( &c[0] );
  };
  void setCoeffs( T* c ) {

    //if ( c.size() != CUBIC_COEFF_ALL ) return;

    for ( int i = 3; i < CUBIC_COEFF_ALL-1; ++i )
      setCoeff( i-3, c[i] );
  };
  void getCoeffs( std::vector<T>& c ) {

//     if ( c.size() != CUBIC_COEFF_ALL ) return;

    for ( int i = 0; i < CUBIC_COEFF; ++i )
      c[i+3] = coeff(i);
  };

private:

  //  0    1    2    3    4    5    6    7    8    9
  // XXX, YYY, ZZZ, XXY, YYZ, ZZX, XYY, YZZ, ZXX, XYZ, 
  //
  // 10  11  12  13  14  15  16  17  18  19
  // XX, YY, ZZ, XY, YZ, ZX, X,  Y,  Z,  const
  //
//   std::vector<T> c_;

#if 0
  T cXXX_, cYYY_, cZZZ_;
  T cXXY_, cYYZ_, cZZX_, cXYY_, cYZZ_, cZXX_, cXYZ_;
  T cXX_, cYY_, cZZ_;
  T cXY_, cYZ_, cZX_;
  T cX_, cY_, cZ_, c0_;
#endif

  T poly_( T x, T y, T z ) const {
    return ( ( (coeff(0) * x + coeff(3) * y + coeff(10)) * x + coeff(6) * y * y + coeff(13) * y + coeff(16) ) * x +
	     ( (coeff(1) * y + coeff(4) * z + coeff(11)) * y + coeff(7) * z * z + coeff(14) * z + coeff(17) ) * y +
	     ( (coeff(2) * z + coeff(5) * x + coeff(12)) * z + coeff(8) * x * x + coeff(15) * x + coeff(18) ) * z + 
	     coeff(9) * x * y * z + coeff(19) );
  };
  
  void polyG_( T g[3], T x, T y, T z ) {
    g[0] = ( (3.0f * coeff(0) * x + 2.0f * (coeff(3) * y + coeff(10))) * x + (coeff(9) * z + coeff(6) * y + coeff(13)) * y + (2.0f * coeff(8) * x + coeff(5) * z + coeff(15)) * z + coeff(16) );
    g[1] = ( (3.0f * coeff(1) * y + 2.0f * (coeff(4) * z + coeff(11))) * y + (coeff(9) * x + coeff(7) * z + coeff(14)) * z + (2.0f * coeff(6) * y + coeff(3) * x + coeff(13)) * x + coeff(17) );
    g[2] = ( (3.0f * coeff(2) * z + 2.0f * (coeff(5) * x + coeff(12))) * z + (coeff(9) * y + coeff(8) * x + coeff(15)) * x + (2.0f * coeff(7) * z + coeff(4) * y + coeff(14)) * y + coeff(18) );
  };

};

typedef Cubic<double> Cubicd;
typedef Cubic<float>  Cubicf;

template <typename T>
class Quadratic : public BasisFunction<T> {

public:

  Quadratic(){
    resize( QUADRATIC_COEFF );
  };

  ~Quadratic(){
    clear();
  };

//   std::vector<T>& coeff() { return c_; };

  void setCoeffs( std::vector<T>& c ) {

//     if ( c.size() != QUADRATIC_COEFF_ALL ) return;
//     for ( int i = 3; i < QUADRATIC_COEFF_ALL-1; ++i )
//       setCoeff( i-3, c[i] );

    setCoeffs( &c[0] );

  };

  void  setCoeffs( T* c ) {

//     if ( c.size() != QUADRATIC_COEFF_ALL ) return;
    
    for ( int i = 3; i < QUADRATIC_COEFF_ALL-1; ++i )
      setCoeff( i-3, c[i] );
  };
  
  void getCoeffs( std::vector<T>& c ) {

    //if ( c.size() != QUADRATIC_COEFF_ALL ) return;

    for ( int i = 0; i < QUADRATIC_COEFF; ++i )
      c[i+3] = coeff(i);
  };

  friend std::ostream& operator <<( std::ostream& o, const Quadratic<T>& q ) {
    return o << "( " << q.coeff(0) << ", " << q.coeff(1) << ", " << q.coeff(2) << ", " << q.coeff(3) << ", " << q.coeff(4) << " " << q.coeff(5) << " " << q.coeff(6) << " " << q.coeff(7) << " " << q.coeff(8) << " " << q.coeff(9) << " )\n";
  };


private:

//   std::vector<T> cxx_[4]; // cXX_, cYY_, cZZ_, .0f
//   std::vector<T> cxy_[4]; // cXY_, cYZ_, cZX_, .0f
//   std::vector<T> cx_[4];  // cX_, cY_, cZ_, c0_

  //
  // 0   1   2   3   4   5   6   7   8   9
  // XX, YY, ZZ, XY, YZ, ZX, X,  Y,  Z,  const
  //
//   std::vector<T> c_;

#if 0
  T cXX_, cYY_, cZZ_;
  T cXY_, cYZ_, cZX_;
  T cX_, cY_, cZ_, c0_;
#endif
  
  T poly_( T x, T y, T z ) const {
    return ( (coeff(0) * x + coeff(3) * y + coeff(6)) * x + 
	     (coeff(1) * y + coeff(4) * z + coeff(7)) * y + 
	     (coeff(2) * z + coeff(5) * x + coeff(8)) * z + coeff(9) );
  };
  
  void polyG_( T g[3], T x, T y, T z ) {
    g[0] = 2.0f * coeff(0) * x + coeff(3) * y + coeff(5) * z + coeff(6);
    g[1] = 2.0f * coeff(1) * y + coeff(4) * z + coeff(3) * x + coeff(7);
    g[2] = 2.0f * coeff(2) * z + coeff(5) * x + coeff(4) * y + coeff(8);
  };

};

typedef Quadratic<double> Quadraticd;
typedef Quadratic<float>  Quadraticf;

template <typename T>
class Linear : public BasisFunction<T> {

public:

  Linear(){
    resize( LINEAR_COEFF );
  };

  ~Linear(){
    clear();
  };

//   std::vector<T>& coeff() { return c_; };

  void setCoeffs( std::vector<T>& c ) {
    
//     if ( c.size() != LINEAR_COEFF_ALL ) return;

    for ( int i = 3; i < LINEAR_COEFF_ALL-1; ++i )
      setCoeff( i-3, c[i] );
  };
  void getCoeffs( std::vector<T>& c ) {

//     if ( c.size() != LINEAR_COEFF_ALL ) return;

    for ( int i = 0; i < LINEAR_COEFF; ++i )
      c[i+3] = coeff(i);
  };

private:

  //
  // 0  1  2  3
  // X, Y, Z, const
  //

//   std::vector<T> c_;

#if 0
  T cX_, cY_, cZ_, c0_;
#endif
  
  T poly_( T x, T y, T z ) const {
    return coeff(0) * x + coeff(1) * y + coeff(2) * z + coeff(3);
  };
  
  void polyG_( T g[3], T x, T y, T z ) {
    g[0] = coeff(0);
    g[1] = coeff(1);
    g[2] = coeff(2);
  };

};

typedef Linear<double> Lineard;
typedef Linear<float>  Linearf;


#endif // _BASISFUNCTION_HXX
