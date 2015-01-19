////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#if defined(WIN32)
#include <float.h>
#endif

#include <vector>
#include <cmath>

#include "Point2.h"
#include "Vector2.h"
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#include "SVD.h"
#include "PointFit2D.hxx"

void PointFit2D::addPointFitElementsDisA( Point2d& p )
{
  double x = p.x;
  double y = p.y;

  addElementDis( 0, 0, 
		 std::pow(x,4) );
  addElementDis( 0, 1,
		 std::pow(x,2)*std::pow(y,2) );
  addElementDis( 0, 2,
		 std::pow(x,3)*y );
  addElementDis( 0, 3,
		 std::pow(x,3) );
  addElementDis( 0, 4,
		 std::pow(x,2)*y );
  addElementDis( 0, 5,
		 std::pow(x,2) );

  addElementDis( 1, 0,
		 std::pow(x,2)*std::pow(y,2) );
  addElementDis( 1, 1,
		 std::pow(y,4) );
  addElementDis( 1, 2, 
		 x*std::pow(y,3) );
  addElementDis( 1, 3, 
		 x*std::pow(y,2) );
  addElementDis( 1, 4, 
		 std::pow(y,3) );
  addElementDis( 1, 5, 
		 std::pow(y,2) );

  addElementDis( 2, 0,  
		 std::pow(x,3)*y );
  addElementDis( 2, 1,  
		 x*std::pow(y,3) );
  addElementDis( 2, 2,  
		 std::pow(x,2)*std::pow(y,2) );
  addElementDis( 2, 3,
		 std::pow(x,2)*y );
  addElementDis( 2, 4,
		 x*std::pow(y,2) );
  addElementDis( 2, 5,
		 x*y );

  addElementDis( 3, 0,
		 std::pow(x,3) );
  addElementDis( 3, 1,
		 x*std::pow(y,2) );
  addElementDis( 3, 2,
		 std::pow(x,2)*y );
  addElementDis( 3, 3,
		 std::pow(x,2) );
  addElementDis( 3, 4,
		 x*y );
  addElementDis( 3, 5,
		 x );

  addElementDis( 4, 0,
		 std::pow(x,2)*y );
  addElementDis( 4, 1,
		 std::pow(y,3) );
  addElementDis( 4, 2,
		 x*std::pow(y,2) );
  addElementDis( 4, 3,
		 x*y );
  addElementDis( 4, 4,
		 std::pow(y,2) );
  addElementDis( 4, 5,
		 y );

  addElementDis( 5, 0,
		 std::pow(x,2) );
  addElementDis( 5, 1,
		 std::pow(y,2) );
  addElementDis( 5, 2,
		 x*y );
  addElementDis( 5, 3,
		 x );
  addElementDis( 5, 4,
		 y );
  addElementDis( 5, 5,
		 1 );
}

void PointFit2D::addPointFitElementsNormA( Point2d& p )
{
  double x = p.x;
  double y = p.y;

  addElementNormA( 0, 0,
		  4*std::pow(x,2) );
//   addElementNormA( 0, 1,
// 		  0 );
  addElementNormA( 0, 2,
		  2*x*y );
  addElementNormA( 0, 3,
		  2*x );
//   addElementNormA( 0, 4,
// 		  0 );
//   addElementNormA( 0, 5,
// 		  0 );
//   addElementNormA( 1, 0,
// 		  0 );
  addElementNormA( 1, 1,
		  4*std::pow(y,2) );
  addElementNormA( 1, 2,
		  2*x*y );
//   addElementNormA( 1, 3,
// 		  0 );
  addElementNormA( 1, 4,
		  2*y );
//   addElementNormA( 1, 5,
// 		  0 );
  addElementNormA( 2, 0,
		  2*x*y );
  addElementNormA( 2, 1,
		  2*x*y );
  addElementNormA( 2, 2,
		  std::pow(x,2) + std::pow(y,2) );
  addElementNormA( 2, 3,
		  y );
  addElementNormA( 2, 4,
		  x );
//   addElementNormA( 2, 5,
// 		  0 );
  addElementNormA( 3, 0,
		  2*x );
//   addElementNormA( 3, 1,
// 		  0 );
  addElementNormA( 3, 2,
		  y );
  addElementNormA( 3, 3,
		  1 );
//   addElementNormA( 3, 4,
// 		  0 );
//   addElementNormA( 3, 5,
// 		  0 );
//   addElementNormA( 4, 0,
// 		  0 );
  addElementNormA( 4, 1,
		  2*y );
  addElementNormA( 4, 2,
		  x );
//   addElementNormA( 4, 3,
// 		  0 );
  addElementNormA( 4, 4,
		  1 );
//   addElementNormA( 4, 5,
// 		  0 );
//   addElementNormA( 5, 0,
// 		  0 );
//   addElementNormA( 5, 1,
// 		  0 );
//   addElementNormA( 5, 2,
// 		  0 );
//   addElementNormA( 5, 3,
// 		  0 );
//   addElementNormA( 5, 4,
// 		  0 );
//   addElementNormA( 5, 5,
// 		  0 );

}

void PointFit2D::addPointFitElementsNormB( Point2d& p, Vector2d& nm )
{
  double x = p.x;
  double y = p.y;
  double nx = nm.x;
  double ny = nm.y;

  addElementNormb( 0, 2*nx*x );
  addElementNormb( 1, 2*ny*y );
  addElementNormb( 2, ny*x + nx*y );
  addElementNormb( 3, nx );
  addElementNormb( 4, ny );
  addElementNormb( 5, 0 );
}

void PointFit2D::addPointFitElementsNormC( Vector2d& nm )
{
  addElementNormc( nm.x * nm.x + nm.y * nm.y );
}

double PointFit2D::error( std::vector<double>& coeff )
{
  if ( !A_ ) return -1;

  double error = .0;

  std::vector<double> mid( NUM_VEC - 1 );
  for ( int i = 1; i < NUM_VEC; ++i ) mid[i-1] = .0;
  for ( int i = 1; i < NUM_VEC; ++i )
    {
      for ( int j = 1; j < NUM_VEC; ++j )
	{
	  mid[i-1] += (A_[i][j] * coeff[j-1]);
	}
    }

  for ( int i = 1; i < NUM_VEC; ++i )
    error += mid[i-1] * coeff[i-1];
  
  for ( int i = 1; i < NUM_VEC; ++i )
    error -= (2 * b_[i] * coeff[i-1]);

  error += c_;
  
  return error;
}

bool PointFit2D::optimize( std::vector<double>& opt )
{
  if ( !A_ ) return false;

  std::vector<double> w( NUM_VEC );
  double** v = new double*[ NUM_VEC ];
  for ( int i = 1; i < NUM_VEC; ++i ) v[i] = new double[ NUM_VEC ];

  SVD::svdcmp( A_, NUM_VEC - 1, NUM_VEC - 1, &w[0], v );

  double wmax=0.0f;
  for (int k = 1; k < NUM_VEC; ++k )
    {
      // 計算に失敗
#if defined(WIN32)
      if ( _isnan( w[k] ) ) return false;
#else
      if ( isnan( w[k] ) ) return false;
#endif
      if ( std::fabs(w[k]) > wmax ) wmax = std::fabs(w[k]);
    }

  double wmin=wmax*0.000001f;
  for (int k = 1; k < NUM_VEC; ++k )
    if ( std::fabs(w[k]) < wmin ) w[k] = 0.;
  
  std::vector<double> x( NUM_VEC );
  SVD::svbksb( A_, &w[0], v, NUM_VEC - 1, NUM_VEC - 1, &b_[0], &x[0] );

  for ( int i = 0; i < NUM_VEC - 1; ++i )
    opt[i] = x[i+1];

  for ( int i = 1; i < NUM_VEC; ++i ) delete v[i];
  delete v;

  return true;
}

// /////////////////////////////////////////////////////////////////////////////////////////////////////

// //
// // PointFit2D から A, b, c を作成
// //
// void PointFit2D::copyPointFitToMatrixVec()
// {
//   int n = 0;
//   for ( int i = 1; i < NUM_VEC; ++i )
//     {
//       for ( int j = i; j < NUM_VEC; ++j )
// 	{
// 	  double t = elmDis_[n] + length_ * elmNorm_[n];
// 	  A_[i][j] = t;
// 	  if ( i != j ) 
// 	    {
// 	      A_[j][i] = t;
// 	    }
// 	  ++n;
// 	}
//     }
//   for ( int i = 1; i < NUM_VEC; ++i )
//     {
//       b_[i] = length_ * elmNorm_[n];
//       ++n;
//     }
  
//   c_ = length_ * elmNorm_[n];
// }

/////////////////////////////////////////////////////////////////////////////////////////////////////
