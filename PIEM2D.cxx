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
#include "PIEM2D.hxx"

void PIEM2D::addPIEMElementsDisA( Point2d& p0, Point2d& p1 )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double x1 = p1.x;
  double y1 = p1.y;

  // A[0][0]
  addElementDis( 0, 
		 (std::pow(x0,4) + std::pow(x0,3)*x1 + std::pow(x0,2)*std::pow(x1,2) + x0*std::pow(x1,3) + 
		  std::pow(x1,4))/5., weightDis_ );

  // A[0][1]
  addElementDis( 1, 
		 (std::pow(x0,2)*(6*std::pow(y0,2) + 3*y0*y1 + std::pow(y1,2)) + 
		  x0*x1*(3*std::pow(y0,2) + 4*y0*y1 + 3*std::pow(y1,2)) + 
		  std::pow(x1,2)*(std::pow(y0,2) + 3*y0*y1 + 6*std::pow(y1,2)))/30., weightDis_ );

  // A[0][2]
  addElementDis( 2, 
		 (std::pow(x0,3)*(4*y0 + y1) + std::pow(x0,2)*x1*(3*y0 + 2*y1) + 
		  x0*std::pow(x1,2)*(2*y0 + 3*y1) + std::pow(x1,3)*(y0 + 4*y1))/20., weightDis_ );

  // A[0][3]
  addElementDis( 3, 
		 (std::pow(x0,3) + std::pow(x0,2)*x1 + x0*std::pow(x1,2) + std::pow(x1,3))/4., weightDis_ );

  // A[0][4]
  addElementDis( 4, 
		 (2*x0*x1*(y0 + y1) + std::pow(x0,2)*(3*y0 + y1) + std::pow(x1,2)*(y0 + 3*y1))/12., weightDis_ );

  // A[0][5]
  addElementDis( 5, 
		 (std::pow(x0,2) + x0*x1 + std::pow(x1,2))/3., weightDis_ );

  // A[1][1]
  addElementDis( 6, 
		 (std::pow(y0,4) + std::pow(y0,3)*y1 + std::pow(y0,2)*std::pow(y1,2) + y0*std::pow(y1,3) + 
		  std::pow(y1,4))/5., weightDis_ );

  // A[1][2]
  addElementDis( 7, 
		 (x0*(4*std::pow(y0,3) + 3*std::pow(y0,2)*y1 + 2*y0*std::pow(y1,2) + std::pow(y1,3)) + 
		  x1*(std::pow(y0,3) + 2*std::pow(y0,2)*y1 + 3*y0*std::pow(y1,2) + 4*std::pow(y1,3)))/
		 20., weightDis_ );

  // A[1][3]
  addElementDis( 8, 
		 (x0*(3*std::pow(y0,2) + 2*y0*y1 + std::pow(y1,2)) + 
		  x1*(std::pow(y0,2) + 2*y0*y1 + 3*std::pow(y1,2)))/12., weightDis_ );

  // A[1][4]
  addElementDis( 9, 
		 (std::pow(y0,3) + std::pow(y0,2)*y1 + y0*std::pow(y1,2) + std::pow(y1,3))/4., weightDis_ );

  // A[1][5]
  addElementDis( 10, 
		 (std::pow(y0,2) + y0*y1 + std::pow(y1,2))/3., weightDis_ );

  // A[2][2]
  addElementDis( 11, 
		 (std::pow(x0,2)*(6*std::pow(y0,2) + 3*y0*y1 + std::pow(y1,2)) + 
		  x0*x1*(3*std::pow(y0,2) + 4*y0*y1 + 3*std::pow(y1,2)) + 
		  std::pow(x1,2)*(std::pow(y0,2) + 3*y0*y1 + 6*std::pow(y1,2)))/30., weightDis_ );

  // A[2][3]
  addElementDis( 12, 
		 (2*x0*x1*(y0 + y1) + std::pow(x0,2)*(3*y0 + y1) + std::pow(x1,2)*(y0 + 3*y1))/12., weightDis_ );

  // A[2][4]
  addElementDis( 13, 
		 (x0*(3*std::pow(y0,2) + 2*y0*y1 + std::pow(y1,2)) + 
		  x1*(std::pow(y0,2) + 2*y0*y1 + 3*std::pow(y1,2)))/12., weightDis_ );

  // A[2][5]
  addElementDis( 14, 
		 (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/6., weightDis_ );

  // A[3][3]
  addElementDis( 15, 
		 (std::pow(x0,2) + x0*x1 + std::pow(x1,2))/3., weightDis_ );

  // A[3][4]
  addElementDis( 16, 
		 (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/6., weightDis_ );

  // A[3][5]
  addElementDis( 17, 
		 (x0 + x1)/2., weightDis_ );

  // A[4][4]
  addElementDis( 18, 
		 (std::pow(y0,2) + y0*y1 + std::pow(y1,2))/3., weightDis_ );

  // A[4][5]
  addElementDis( 19, 
		 (y0 + y1)/2., weightDis_ );

  // A[5][5]
  addElementDis( 20, 
		 1, weightDis_ );

//   // A[1][0]
//   addElementDis( 0, 
// (std::pow(x0,2)*(6*std::pow(y0,2) + 3*y0*y1 + std::pow(y1,2)) + 
//     x0*x1*(3*std::pow(y0,2) + 4*y0*y1 + 3*std::pow(y1,2)) + 
//     std::pow(x1,2)*(std::pow(y0,2) + 3*y0*y1 + 6*std::pow(y1,2)))/30.

//   // A[2][0]
//   addElementDis( 0, 
// (std::pow(x0,3)*(4*y0 + y1) + std::pow(x0,2)*x1*(3*y0 + 2*y1) + 
//     x0*std::pow(x1,2)*(2*y0 + 3*y1) + std::pow(x1,3)*(y0 + 4*y1))/20.

//   // A[2][1]
//   addElementDis( 0, 
// (x0*(4*std::pow(y0,3) + 3*std::pow(y0,2)*y1 + 2*y0*std::pow(y1,2) + std::pow(y1,3)) + 
//     x1*(std::pow(y0,3) + 2*std::pow(y0,2)*y1 + 3*y0*std::pow(y1,2) + 4*std::pow(y1,3)))/
//   20.

//   // A[3][0]
//   addElementDis( 0, 
// (std::pow(x0,3) + std::pow(x0,2)*x1 + x0*std::pow(x1,2) + std::pow(x1,3))/4.

//   // A[3][1]
//   addElementDis( 0, 
// (x0*(3*std::pow(y0,2) + 2*y0*y1 + std::pow(y1,2)) + 
//     x1*(std::pow(y0,2) + 2*y0*y1 + 3*std::pow(y1,2)))/12.

//   // A[3][2]
//   addElementDis( 0, 
// (2*x0*x1*(y0 + y1) + std::pow(x0,2)*(3*y0 + y1) + std::pow(x1,2)*(y0 + 3*y1))/12.

//   // A[4][0]
//   addElementDis( 0, 
// (2*x0*x1*(y0 + y1) + std::pow(x0,2)*(3*y0 + y1) + std::pow(x1,2)*(y0 + 3*y1))/12.

//   // A[4][1]
//   addElementDis( 0, 
// (std::pow(y0,3) + std::pow(y0,2)*y1 + y0*std::pow(y1,2) + std::pow(y1,3))/4.

//   // A[4][2]
//   addElementDis( 0, 
// (x0*(3*std::pow(y0,2) + 2*y0*y1 + std::pow(y1,2)) + 
//     x1*(std::pow(y0,2) + 2*y0*y1 + 3*std::pow(y1,2)))/12.

//   // A[4][3]
//   addElementDis( 0, 
// (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/6.

//   // A[5][0]
//   addElementDis( 0, 
// (std::pow(x0,2) + x0*x1 + std::pow(x1,2))/3.

//   // A[5][1]
//   addElementDis( 0, 
// (std::pow(y0,2) + y0*y1 + std::pow(y1,2))/3.

//   // A[5][2]
//   addElementDis( 0, 
// (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/6.

//   // A[5][3]
//   addElementDis( 0, 
// (x0 + x1)/2.

//   // A[5][4]
//   addElementDis( 0, 
// (y0 + y1)/2.

}

void PIEM2D::addPIEMElementsNormA( Point2d& p0, Point2d& p1 )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double x1 = p1.x;
  double y1 = p1.y;
  
  // 0
  // A[0][0]
  addElementNorm( 0, 
		  (4*(std::pow(x0,2) + x0*x1 + std::pow(x1,2)))/3., weightNrm_ );
  // 1
  // A[0][1]
//   addElementNorm( 1, 
// 		  0, weightNrm_ );
  // 2
  // A[0][2]
  addElementNorm( 2, 
		  (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/3., weightNrm_ );
  // 3
  // A[0][3]
  addElementNorm( 3, 
		  x0 + x1, weightNrm_ );
  // 4
  // A[0][4]
//   addElementNorm( 4, 
// 		  0, weightNrm_ );
  // 5
  // A[0][5]
//   addElementNorm( 5, 
// 		  0, weightNrm_ );

  // 6
  // A[1][1]
  addElementNorm( 6, 
		  (4*(std::pow(y0,2) + y0*y1 + std::pow(y1,2)))/3., weightNrm_ );
  // 7
  // A[1][2]
  addElementNorm( 7, 
		  (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/3., weightNrm_ );
  // 8
  // A[1][3]
//   addElementNorm( 8, 
// 		  0, weightNrm_ );
  // 9
  // A[1][4]
  addElementNorm( 9, 
		  y0 + y1, weightNrm_ );
  // 10
  // A[1][5]
//   addElementNorm( 10, 
// 		  0, weightNrm_ );

  // 11
  // A[2][2]
  addElementNorm( 11,
		  (std::pow(x0,2) + x0*x1 + std::pow(x1,2))/3. + 
		  (std::pow(y0,2) + y0*y1 + std::pow(y1,2))/3., weightNrm_ );
  // 12
  // A[2][3]
  addElementNorm( 12, 
		  (y0 + y1)/2., weightNrm_ );
  // 13
  // A[2][4]
  addElementNorm( 13, 
		  (x0 + x1)/2., weightNrm_ );
  // 14
  // A[2][5]
//   addElementNorm( 14, 
// 		  0, weightNrm_ );

  // 15
  // A[3][3]
  addElementNorm( 15, 
		  1, weightNrm_ );
  // 16
  // A[3][4]
//   addElementNorm( 16, 
// 		  0, weightNrm_ );
  // 17
  // A[3][5]
//   addElementNorm( 17, 
// 		  0, weightNrm_ );

  // 18
  // A[4][4]
  addElementNorm( 18, 
		  1, weightNrm_ );
  // 19
  // A[4][5]
//   addElementNorm( 19, 
// 		  0, weightNrm_ );

  // 20
  // A[5][5]
//   addElementNorm( 20, 
// 		  0, weightNrm_ );

//   // A[1][0]
// 0

//   // A[2][0]
// (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/3.
//   // A[2][1]
// (x0*(2*y0 + y1) + x1*(y0 + 2*y1))/3.

//   // A[3][0]
// x0 + x1
//   // A[3][1]
// 0
//   // A[3][2]
// (y0 + y1)/2.

//   // A[4][0]
// 0
//   // A[4][1]
// y0 + y1
//   // A[4][2]
// (x0 + x1)/2.
//   // A[4][3]
// 0

//   // A[5][0]
// 0
//   // A[5][1]
// 0
//   // A[5][2]
// 0
//   // A[5][3]
// 0
//   // A[5][4]
// 0
}

void PIEM2D::addPIEMElementsNormB( Point2d& p0, Point2d& p1, Vector2d& nm )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double x1 = p1.x;
  double y1 = p1.y;
  double nx = nm.x;
  double ny = nm.y;

  addElementNorm( 21,
		  nx*(x0 + x1), weightNrm_ );
  addElementNorm( 22, 
		  ny*(y0 + y1), weightNrm_ );
  addElementNorm( 23, 
		  (ny*(x0 + x1))/2. + (nx*(y0 + y1))/2., weightNrm_ );
  addElementNorm( 24, 
		  nx, weightNrm_ );
  addElementNorm( 25, 
		  ny, weightNrm_ );
//   addElementNorm( 26, 
// 		  0, weightNrm_ );

}

void PIEM2D::addPIEMElementsNormC( Vector2d& nm )
{
  addElementNorm( 27,
		  nm.x * nm.x + nm.y * nm.y, weightNrm_ );
}

double PIEM2D::error( std::vector<double>& coeff )
{
  if ( !A_ ) initMatrixVec();

  copyPIEMToMatrixVec();

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
  
  clearMatrixVec();

  return error;
}

bool PIEM2D::optimize( std::vector<double>& opt )
{
  if ( !A_ ) initMatrixVec();

  // set Matrix and Vector
  copyPIEMToMatrixVec();

  std::vector<double> w( NUM_VEC );
  double** v = new double*[ NUM_VEC ];
  for ( int i = 1; i < NUM_VEC; ++i ) v[i] = new double[ NUM_VEC ];

//   std::vector< std::vector<double> > v( NUM_VEC );
// //   v.resize( NUM_VEC );
//   for ( int i = 1; i < NUM_VEC; ++i ) v[i].resize( NUM_VEC );

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

  clearMatrixVec();

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

//
// PIEM2D から A, b, c を作成
//
void PIEM2D::copyPIEMToMatrixVec()
{
  int n = 0;
  for ( int i = 1; i < NUM_VEC; ++i )
    {
      for ( int j = i; j < NUM_VEC; ++j )
	{
	  double t = elmDis_[n] + length_ * elmNorm_[n];
	  A_[i][j] = t;
	  if ( i != j ) 
	    {
	      A_[j][i] = t;
	    }
	  ++n;
	}
    }
  for ( int i = 1; i < NUM_VEC; ++i )
    {
      b_[i] = length_ * elmNorm_[n];
      ++n;
    }
  
  c_ = length_ * elmNorm_[n];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////
