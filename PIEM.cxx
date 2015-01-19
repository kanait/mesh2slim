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

#include <cmath>

using namespace std;

#include "NRSVD.h"
#include "PIEM.hxx"

#if _USE_MKL_
#include <itpp/itbase.h>
using namespace itpp;
#endif // _USE_MKL_

void PIEM::addPIEMElementsDisA( Point3d& p0, Point3d& p1, Point3d& p2 )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double z0 = p0.z;
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;
  double x2 = p2.x;
  double y2 = p2.y;
  double z2 = p2.z;

  double p4x0 = x0*x0*x0*x0;
  double p3x0 = x0*x0*x0;
  double p2x0 = x0*x0;
  double p4x1 = x1*x1*x1*x1;
  double p3x1 = x1*x1*x1;
  double p2x1 = x1*x1;
  double p4x2 = x2*x2*x2*x2;
  double p3x2 = x2*x2*x2;
  double p2x2 = x2*x2;

  double p4y0 = y0*y0*y0*y0;
  double p3y0 = y0*y0*y0;
  double p2y0 = y0*y0;
  double p4y1 = y1*y1*y1*y1;
  double p3y1 = y1*y1*y1;
  double p2y1 = y1*y1;
  double p4y2 = y2*y2*y2*y2;
  double p3y2 = y2*y2*y2;
  double p2y2 = y2*y2;

  double p4z0 = z0*z0*z0*z0;
  double p3z0 = z0*z0*z0;
  double p2z0 = z0*z0;
  double p4z1 = z1*z1*z1*z1;
  double p3z1 = z1*z1*z1;
  double p2z1 = z1*z1;
  double p4z2 = z2*z2*z2*z2;
  double p3z2 = z2*z2*z2;
  double p2z2 = z2*z2;

  // 0
  //   A_[0][0] = 
  addElementDis( 0, 
	      (p4x0 + p4x1 + p3x1*x2 + p2x1*p2x2 + 
	       x1*p3x2 + p4x2 + p3x0*(x1 + x2) + 
	       x0*(x1 + x2)*(p2x1 + p2x2) + 
	       p2x0*(p2x1 + x1*x2 + p2x2))/30., weightDis_ );

  // 1
  //   A_[0][1] = A_[1][0] = 
  addElementDis( 1, 
	      (p2x2*(p2y0 + y0*y1 + p2y1 + 3*(y0 + y1)*y2 + 
		     6*p2y2) + x1*x2*
	       (p2y0 + 3*p2y1 + 4*y1*y2 + 3*p2y2 + 
		2*y0*(y1 + y2)) + p2x0*
	       (6*p2y0 + p2y1 + y1*y2 + p2y2 + 3*y0*(y1 + y2)) + 
	       p2x1*(p2y0 + 6*p2y1 + 3*y1*y2 + p2y2 + 
		     y0*(3*y1 + y2)) 
	       + x0*(x1* (3*p2y0 + 4*y0*y1 + 3*p2y1 + 2*(y0 + y1)*y2 + 
			  p2y2) + x2*(3*p2y0 + p2y1 + 2*y1*y2 + 
				      3*p2y2 + 2*y0*(y1 + 2*y2))))/180., weightDis_ );

  // 2
  //   A_[0][2] = A_[2][0] = 
  addElementDis( 2, 
	      (p2x2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		     6*p2z2) + x1*x2*
	       (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + z2)) + p2x0*
	       (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	       p2x1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		     z0*(3*z1 + z2)) 
	       + x0*(x1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			  p2z2) + x2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				      3*p2z2 + 2*z0*(z1 + 2*z2))))/180., weightDis_ );

  // 3
  //   A_[0][3] = A_[3][0] = 
  addElementDis( 3, 
	      (p3x0*(4*y0 + y1 + y2) + p3x1*(y0 + 4*y1 + y2) + 
	       p2x1*x2*(y0 + 3*y1 + 2*y2) + x1*p2x2*(y0 + 2*y1 + 3*y2) + 
	       p3x2*(y0 + y1 + 4*y2) + 
	       p2x0*(x1*(3*y0 + 2*y1 + y2) + x2*(3*y0 + y1 + 2*y2)) + 
	       x0*(2*x1*x2*(y0 + y1 + y2) + p2x1*(2*y0 + 3*y1 + y2) + 
		   p2x2*(2*y0 + y1 + 3*y2)))/120., weightDis_ );

  // 4
  //   A_[0][4] = A_[4][0] = 
  addElementDis( 4, 
	      (p2x0*(3*y0*(4*z0 + z1 + z2) + y1*(3*z0 + 2*z1 + z2) + 
		     y2*(3*z0 + z1 + 2*z2)) + 
	       p2x1*(y0*(2*z0 + 3*z1 + z2) + 3*y1*(z0 + 4*z1 + z2) + 
		     y2*(z0 + 3*z1 + 2*z2)) + 
	       2*x1*x2*(y0*(z0 + z1 + z2) + y1*(z0 + 3*z1 + 2*z2) + 
			y2*(z0 + 2*z1 + 3*z2)) + 
	       p2x2*(y0*(2*z0 + z1 + 3*z2) + y1*(z0 + 2*z1 + 3*z2) + 
		     3*y2*(z0 + z1 + 4*z2)) + 
	       2*x0*(x1*(y2*(z0 + z1 + z2) + y0*(3*z0 + 2*z1 + z2) + 
			 y1*(2*z0 + 3*z1 + z2)) + 
		     x2*(y1*(z0 + z1 + z2) + y0*(3*z0 + z1 + 2*z2) + 
			 y2*(2*z0 + z1 + 3*z2))))/360., weightDis_ );

  // 5
  //   A_[0][5] = A_[5][0] = 
  addElementDis( 5, 
	      (p3x0*(4*z0 + z1 + z2) + p3x1*(z0 + 4*z1 + z2) + 
	       p2x1*x2*(z0 + 3*z1 + 2*z2) + x1*p2x2*(z0 + 2*z1 + 3*z2) + 
	       p3x2*(z0 + z1 + 4*z2) + 
	       p2x0*(x1*(3*z0 + 2*z1 + z2) + x2*(3*z0 + z1 + 2*z2)) + 
	       x0*(2*x1*x2*(z0 + z1 + z2) + p2x1*(2*z0 + 3*z1 + z2) + 
		   p2x2*(2*z0 + z1 + 3*z2)))/120., weightDis_ );

  // 6
  //   A_[0][6] = A_[6][0] = 
  addElementDis( 6, 
	      (p3x0 + p2x0*(x1 + x2) + 
	       (x1 + x2)*(p2x1 + p2x2) + 
	       x0*(p2x1 + x1*x2 + p2x2))/20., weightDis_ );
  
  // 7
  //   A_[0][7] = A_[7][0] = 
  addElementDis( 7, 
	      (p2x0*(3*y0 + y1 + y2) + p2x1*(y0 + 3*y1 + y2) + 
	       p2x2*(y0 + y1 + 3*y2) + x1*x2*(y0 + 2*(y1 + y2)) + 
	       x0*(x1*(2*y0 + 2*y1 + y2) + x2*(2*y0 + y1 + 2*y2)))/60., weightDis_ );

  // 8
  //   A_[0][8] = A_[8][0] = 
  addElementDis( 8, 
	      (p2x0*(3*z0 + z1 + z2) + p2x1*(z0 + 3*z1 + z2) + 
	       p2x2*(z0 + z1 + 3*z2) + x1*x2*(z0 + 2*(z1 + z2)) + 
	       x0*(x1*(2*z0 + 2*z1 + z2) + x2*(2*z0 + z1 + 2*z2)))/60., weightDis_ );

  // 9
  //   A_[0][9] = A_[9][0] = 
  addElementDis( 9, 
	      (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12., weightDis_ );


  // 10
  //   A_[1][1] = 
  addElementDis( 10, 
	      (p4y0 + p4y1 + p3y1*y2 + p2y1*p2y2 + 
	       y1*p3y2 + p4y2 + p3y0*(y1 + y2) + 
	       y0*(y1 + y2)*(p2y1 + p2y2) + 
	       p2y0*(p2y1 + y1*y2 + p2y2))/30., weightDis_ );

  // 11
  //   A_[1][2] = A_[2][1] = 
  addElementDis( 11, 
	      (p2y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		     6*p2z2) + y1*y2*
	       (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + z2)) + p2y0*
	       (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	       p2y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		     z0*(3*z1 + z2)) 
	       + y0*(y1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			  p2z2) + y2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				      3*p2z2 + 2*z0*(z1 + 2*z2))))/180., weightDis_ );

  // 12
  //   A_[1][3] = A_[3][1] = 
  addElementDis( 12, 
	      (x2*(p3y0 + p2y0*y1 + y0*p2y1 + p3y1 + 
		   2*(p2y0 + y0*y1 + p2y1)*y2 + 
		   3*(y0 + y1)*p2y2 + 4*p3y2) + 
	       x0*(4*p3y0 + 3*p2y0*(y1 + y2) + 
		   (y1 + y2)*(p2y1 + p2y2) + 
		   2*y0*(p2y1 + y1*y2 + p2y2)) + 
	       x1*(p3y0 + 4*p3y1 + 3*p2y1*y2 + 2*y1*p2y2 + 
		   p3y2 + p2y0*(2*y1 + y2) + 
		   y0*(3*p2y1 + 2*y1*y2 + p2y2)))/120., weightDis_ );

  // 13
  //   A_[1][4] = A_[4][1] = 
  addElementDis( 13, 
	      (p3y0*(4*z0 + z1 + z2) + p3y1*(z0 + 4*z1 + z2) + 
	       p2y1*y2*(z0 + 3*z1 + 2*z2) + y1*p2y2*(z0 + 2*z1 + 3*z2) + 
	       p3y2*(z0 + z1 + 4*z2) + 
	       p2y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2)) + 
	       y0*(2*y1*y2*(z0 + z1 + z2) + p2y1*(2*z0 + 3*z1 + z2) + 
		   p2y2*(2*z0 + z1 + 3*z2)))/120., weightDis_ );

  // 14
  //   A_[1][5] = A_[5][1] = 
  addElementDis( 14, 
	      (x1*(p2y0*(3*z0 + 2*z1 + z2) + 3*p2y1*(z0 + 4*z1 + z2) + 
		   2*y1*y2*(z0 + 3*z1 + 2*z2) + p2y2*(z0 + 2*z1 + 3*z2) + 
		   2*y0*(y2*(z0 + z1 + z2) + y1*(2*z0 + 3*z1 + z2))) + 
	       x0*(2*y1*y2*(z0 + z1 + z2) + 3*p2y0*(4*z0 + z1 + z2) + 
		   p2y1*(2*z0 + 3*z1 + z2) + p2y2*(2*z0 + z1 + 3*z2) + 
		   2*y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2))) + 
	       x2*(p2y0*(3*z0 + z1 + 2*z2) + p2y1*(z0 + 3*z1 + 2*z2) + 
		   2*y1*y2*(z0 + 2*z1 + 3*z2) + 3*p2y2*(z0 + z1 + 4*z2) + 
		   2*y0*(y1*(z0 + z1 + z2) + y2*(2*z0 + z1 + 3*z2))))/360., weightDis_ );

  // 15
  //   A_[1][6] = A_[6][1] = 
  addElementDis( 15, 
	      (x1*(p2y0 + 2*y0*y1 + 3*p2y1 + y0*y2 + 2*y1*y2 + 
		   p2y2) + x2*(p2y0 + y0*y1 + p2y1 + 2*y0*y2 + 
			       2*y1*y2 + 3*p2y2) + 
	       x0*(3*p2y0 + p2y1 + y1*y2 + p2y2 + 2*y0*(y1 + y2)))/60., weightDis_ );

  // 16
  //   A_[1][7] = A_[7][1] = 
  addElementDis( 16, 
	      (p3y0 + p2y0*(y1 + y2) + 
	       (y1 + y2)*(p2y1 + p2y2) + 
	       y0*(p2y1 + y1*y2 + p2y2))/20., weightDis_ );

  // 17
  //   A_[1][8] = A_[8][1] = 
  addElementDis( 17, 
	      (p2y0*(3*z0 + z1 + z2) + p2y1*(z0 + 3*z1 + z2) + 
	       p2y2*(z0 + z1 + 3*z2) + y1*y2*(z0 + 2*(z1 + z2)) + 
	       y0*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2)))/60., weightDis_ );

  // 18
  //   A_[1][9] = A_[9][1] = 
  addElementDis( 18, 
	      (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12., weightDis_ );


  // 19
  //   A_[2][2] = 
  addElementDis( 19, 
	      (p4z0 + p4z1 + p3z1*z2 + p2z1*p2z2 + 
	       z1*p3z2 + p4z2 + p3z0*(z1 + z2) + 
	       z0*(z1 + z2)*(p2z1 + p2z2) + 
	       p2z0*(p2z1 + z1*z2 + p2z2))/30., weightDis_ );

  // 20
  //   A_[2][3] = A_[3][2] = 
  addElementDis( 20, 
	      (x1*(y0*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		       p2z2) + y2*(p2z0 + 3*p2z1 + 4*z1*z2 + 
				   3*p2z2 + 2*z0*(z1 + z2)) + 
		   2*y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
			 z0*(3*z1 + z2))) + x2*
	       (2*y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		      6*p2z2) + y1*
		(p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		 2*z0*(z1 + z2)) + y0*
		(3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		 2*z0*(z1 + 2*z2))) + 
	       x0*(y1*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		       p2z2) + 2*y0*
		   (6*p2z0 + p2z1 + z1*z2 + p2z2 + 
		    3*z0*(z1 + z2)) + y2*
		   (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		    2*z0*(z1 + 2*z2))))/360., weightDis_ );

  // 21
  //   A_[2][4] = A_[4][2] = 
  addElementDis( 21, 
	      (y2*(p3z0 + p2z0*z1 + z0*p2z1 + p3z1 + 
		   2*(p2z0 + z0*z1 + p2z1)*z2 + 
		   3*(z0 + z1)*p2z2 + 4*p3z2) + 
	       y0*(4*p3z0 + 3*p2z0*(z1 + z2) + 
		   (z1 + z2)*(p2z1 + p2z2) + 
		   2*z0*(p2z1 + z1*z2 + p2z2)) + 
	       y1*(p3z0 + 4*p3z1 + 3*p2z1*z2 + 2*z1*p2z2 + 
		   p3z2 + p2z0*(2*z1 + z2) + 
		   z0*(3*p2z1 + 2*z1*z2 + p2z2)))/120., weightDis_ );

  // 22
  //   A_[2][5] = A_[5][2] = 
  addElementDis( 22, 
	      (x2*(p3z0 + p2z0*z1 + z0*p2z1 + p3z1 + 
		   2*(p2z0 + z0*z1 + p2z1)*z2 + 
		   3*(z0 + z1)*p2z2 + 4*p3z2) + 
	       x0*(4*p3z0 + 3*p2z0*(z1 + z2) + 
		   (z1 + z2)*(p2z1 + p2z2) + 
		   2*z0*(p2z1 + z1*z2 + p2z2)) + 
	       x1*(p3z0 + 4*p3z1 + 3*p2z1*z2 + 2*z1*p2z2 + 
		   p3z2 + p2z0*(2*z1 + z2) + 
		   z0*(3*p2z1 + 2*z1*z2 + p2z2)))/120., weightDis_ );

  // 23
  //   A_[2][6] = A_[6][2] = 
  addElementDis( 23, 
	      (x1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		   p2z2) + x2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			       2*z1*z2 + 3*p2z2) + 
	       x0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60., weightDis_ );

  // 24
  //   A_[2][7] = A_[7][2] = 
  addElementDis( 24, 
	      (y1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		   p2z2) + y2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			       2*z1*z2 + 3*p2z2) + 
	       y0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60., weightDis_ );

  // 25
  //   A_[2][8] = A_[8][2] = 
  addElementDis( 25, 
	      (p3z0 + p2z0*(z1 + z2) + 
	       (z1 + z2)*(p2z1 + p2z2) + 
	       z0*(p2z1 + z1*z2 + p2z2))/20., weightDis_ );

  // 26
  //   A_[2][9] = A_[9][2] = 
  addElementDis( 26, 
	      (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12., weightDis_ );


  // 27
  //   A_[3][3] = 
  addElementDis( 27, 
	      (p2x2*(p2y0 + y0*y1 + p2y1 + 3*(y0 + y1)*y2 + 
		     6*p2y2) + x1*x2*
	       (p2y0 + 3*p2y1 + 4*y1*y2 + 3*p2y2 + 
		2*y0*(y1 + y2)) + p2x0*
	       (6*p2y0 + p2y1 + y1*y2 + p2y2 + 3*y0*(y1 + y2)) + 
	       p2x1*(p2y0 + 6*p2y1 + 3*y1*y2 + p2y2 + 
		     y0*(3*y1 + y2)) 
	       + x0*(x1* (3*p2y0 + 4*y0*y1 + 3*p2y1 + 2*(y0 + y1)*y2 + 
			  p2y2) + x2*(3*p2y0 + p2y1 + 2*y1*y2 + 
				      3*p2y2 + 2*y0*(y1 + 2*y2))))/180., weightDis_ );

  // 28
  //   A_[3][4] = A_[4][3] = 
  addElementDis( 28, 
	      (x1*(p2y0*(3*z0 + 2*z1 + z2) + 3*p2y1*(z0 + 4*z1 + z2) + 
		   2*y1*y2*(z0 + 3*z1 + 2*z2) + p2y2*(z0 + 2*z1 + 3*z2) + 
		   2*y0*(y2*(z0 + z1 + z2) + y1*(2*z0 + 3*z1 + z2))) + 
	       x0*(2*y1*y2*(z0 + z1 + z2) + 3*p2y0*(4*z0 + z1 + z2) + 
		   p2y1*(2*z0 + 3*z1 + z2) + p2y2*(2*z0 + z1 + 3*z2) + 
		   2*y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2))) + 
	       x2*(p2y0*(3*z0 + z1 + 2*z2) + p2y1*(z0 + 3*z1 + 2*z2) + 
		   2*y1*y2*(z0 + 2*z1 + 3*z2) + 3*p2y2*(z0 + z1 + 4*z2) + 
		   2*y0*(y1*(z0 + z1 + z2) + y2*(2*z0 + z1 + 3*z2))))/360., weightDis_ );

  // 29
  //   A_[3][5] = A_[5][3] = 
  addElementDis( 29, 
	      (p2x0*(3*y0*(4*z0 + z1 + z2) + y1*(3*z0 + 2*z1 + z2) + 
		     y2*(3*z0 + z1 + 2*z2)) + 
	       p2x1*(y0*(2*z0 + 3*z1 + z2) + 3*y1*(z0 + 4*z1 + z2) + 
		     y2*(z0 + 3*z1 + 2*z2)) + 
	       2*x1*x2*(y0*(z0 + z1 + z2) + y1*(z0 + 3*z1 + 2*z2) + 
			y2*(z0 + 2*z1 + 3*z2)) + 
	       p2x2*(y0*(2*z0 + z1 + 3*z2) + y1*(z0 + 2*z1 + 3*z2) + 
		     3*y2*(z0 + z1 + 4*z2)) + 
	       2*x0*(x1*(y2*(z0 + z1 + z2) + y0*(3*z0 + 2*z1 + z2) + 
			 y1*(2*z0 + 3*z1 + z2)) + 
		     x2*(y1*(z0 + z1 + z2) + y0*(3*z0 + z1 + 2*z2) + 
			 y2*(2*z0 + z1 + 3*z2))))/360., weightDis_ );

  // 30
  //   A_[3][6] = A_[6][3] = 
  addElementDis( 30, 
	      (p2x0*(3*y0 + y1 + y2) + p2x1*(y0 + 3*y1 + y2) + 
	       p2x2*(y0 + y1 + 3*y2) + x1*x2*(y0 + 2*(y1 + y2)) + 
	       x0*(x1*(2*y0 + 2*y1 + y2) + x2*(2*y0 + y1 + 2*y2)))/60., weightDis_ );

  // 31
  //   A_[3][7] = A_[7][3] = 
  addElementDis( 31, 
	      (x1*(p2y0 + 2*y0*y1 + 3*p2y1 + y0*y2 + 2*y1*y2 + 
		   p2y2) + x2*(p2y0 + y0*y1 + p2y1 + 2*y0*y2 + 
			       2*y1*y2 + 3*p2y2) + 
	       x0*(3*p2y0 + p2y1 + y1*y2 + p2y2 + 2*y0*(y1 + y2)))/60., weightDis_ );

  // 32
  //   A_[3][8] = A_[8][3] = 
  addElementDis( 32, 
	      (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
		   y2*(2*z0 + z1 + 2*z2)) + 
	       x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
		   y1*(z0 + 2*(z1 + z2))) + 
	       x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
		   y2*(z0 + 2*(z1 + z2))))/120., weightDis_ );

  // 33
  //   A_[3][9] = A_[9][3] = 
  addElementDis( 33, 
	      (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24., weightDis_ );


  // 34
  //   A_[4][4] = 
  addElementDis( 34, 
	      (p2y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		     6*p2z2) + y1*y2*
	       (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + z2)) + p2y0*
	       (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	       p2y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		     z0*(3*z1 + z2)) 
	       + y0*(y1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			  p2z2) + y2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				      3*p2z2 + 2*z0*(z1 + 2*z2))))/180., weightDis_ );

  // 35
  //   A_[4][5] = A_[5][4] = 
  addElementDis( 35, 
	      (x1*(y0*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		       p2z2) + y2*(p2z0 + 3*p2z1 + 4*z1*z2 + 
				   3*p2z2 + 2*z0*(z1 + z2)) + 
		   2*y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
			 z0*(3*z1 + z2))) + x2*
	       (2*y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		      6*p2z2) + y1*
		(p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		 2*z0*(z1 + z2)) + y0*
		(3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		 2*z0*(z1 + 2*z2))) + 
	       x0*(y1*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		       p2z2) + 2*y0*
		   (6*p2z0 + p2z1 + z1*z2 + p2z2 + 
		    3*z0*(z1 + z2)) + y2*
		   (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		    2*z0*(z1 + 2*z2))))/360., weightDis_ );

  // 36
  //   A_[4][6] = A_[6][4] = 
  addElementDis( 36, 
	      (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
		   y2*(2*z0 + z1 + 2*z2)) + 
	       x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
		   y1*(z0 + 2*(z1 + z2))) + 
	       x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
		   y2*(z0 + 2*(z1 + z2))))/120., weightDis_ );

  // 37
  //   A_[4][7] = A_[7][4] = 
  addElementDis( 37, 
	      (p2y0*(3*z0 + z1 + z2) + p2y1*(z0 + 3*z1 + z2) + 
	       p2y2*(z0 + z1 + 3*z2) + y1*y2*(z0 + 2*(z1 + z2)) + 
	       y0*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2)))/60., weightDis_ );

  // 38
  //   A_[4][8] = A_[8][4] = 
  addElementDis( 38, 
	      (y1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		   p2z2) + y2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			       2*z1*z2 + 3*p2z2) + 
	       y0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60., weightDis_ );

  // 39
  //   A_[4][9] = A_[9][4] = 
  addElementDis( 39, 
	      (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24., weightDis_ );


  // 40
  //   A_[5][5] = 
  addElementDis( 40, 
	      (p2x2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		     6*p2z2) + x1*x2*
	       (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + z2)) + p2x0*
	       (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	       p2x1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		     z0*(3*z1 + z2)) 
	       + x0*(x1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			  p2z2) + x2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				      3*p2z2 + 2*z0*(z1 + 2*z2))))/180., weightDis_ );

  // 41
  //   A_[5][6] = A_[6][5] = 
  addElementDis( 41, 
	      (p2x0*(3*z0 + z1 + z2) + p2x1*(z0 + 3*z1 + z2) + 
	       p2x2*(z0 + z1 + 3*z2) + x1*x2*(z0 + 2*(z1 + z2)) + 
	       x0*(x1*(2*z0 + 2*z1 + z2) + x2*(2*z0 + z1 + 2*z2)))/60., weightDis_ );

  // 42
  //   A_[5][7] = A_[7][5] = 
  addElementDis( 42, 
	      (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
		   y2*(2*z0 + z1 + 2*z2)) + 
	       x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
		   y1*(z0 + 2*(z1 + z2))) + 
	       x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
		   y2*(z0 + 2*(z1 + z2))))/120., weightDis_ );

  // 43
  //   A_[5][8] = A_[8][5] = 
  addElementDis( 43, 
	      (x1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		   p2z2) + x2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			       2*z1*z2 + 3*p2z2) + 
	       x0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60., weightDis_ );

  // 44
  //   A_[5][9] = A_[9][5] = 
  addElementDis( 44, 
	      (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24., weightDis_ );

  // 45
  //   A_[6][6] = 
  addElementDis( 45, 
	      (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12., weightDis_ );

  // 46
  //   A_[6][7] = A_[7][6] = 
  addElementDis( 46, 
	      (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24., weightDis_ );

  // 47
  //   A_[6][8] = A_[8][6] = 
  addElementDis( 47, 
	      (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24., weightDis_ );

  // 48
  //   A_[6][9] = A_[9][6] = 
  addElementDis( 48, 
	      (x0 + x1 + x2)/6., weightDis_ );

  // 49
  //   A_[7][7] = 
  addElementDis( 49, 
	      (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12., weightDis_ );

  // 50
  //   A_[7][8] = A_[8][7] = 
  addElementDis( 50, 
	      (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24., weightDis_ );

  // 51
  //   A_[7][9] = A_[9][7] = 
  addElementDis( 51, 
	      (y0 + y1 + y2)/6., weightDis_ );

  // 52
  //   A_[8][8] = 
  addElementDis( 52, 
	      (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12., weightDis_ );

  // 53
  //   A_[8][9] = A_[9][8] = 
  addElementDis( 53, 
	      (z0 + z1 + z2)/6., weightDis_ );

  // 54
  //   A_[9][9] = 
  addElementDis( 54, 
	      0.5, weightDis_ );
}

void PIEM::addPIEMElementsNormA( Point3d& p0, Point3d& p1, Point3d& p2 )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double z0 = p0.z;
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;
  double x2 = p2.x;
  double y2 = p2.y;
  double z2 = p2.z;
  
  double p2x0 = x0*x0;
  double p2x1 = x1*x1;
  double p2x2 = x2*x2;

  double p2y0 = y0*y0;
  double p2y1 = y1*y1;
  double p2y2 = y2*y2;

  double p2z0 = z0*z0;
  double p2z1 = z1*z1;
  double p2z2 = z2*z2;

//   double t;

  // 0
  //   A_[0][0]+=
  addElementNorm( 0, 
	      (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/3., weightNrm_ );

  // 1
  //   A_[0][1]+=0;

  // 2
  //   A_[0][2]+=0;

  // 3
  //   A_[0][3]+=t;
  //   A_[3][0]+=t;
  addElementNorm( 3, 
	      (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12., weightNrm_ );

  // 4
  //   A_[0][4]+=0;

  // 5
  //   A_[0][5]+=t;
  //   A_[5][0]+=t;
  addElementNorm( 5, 
	      (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12., weightNrm_ );

  // 6
  //   A_[0][6]+=t;
  //   A_[6][0]+=t;
  addElementNorm( 6, 
	      (x0 + x1 + x2)/3., weightNrm_ );

  // 7
  //   A_[0][7]+=0;

  // 8
  //   A_[0][8]+=0;

  // 9
  //   A_[0][9]+=0;

  //   A_[1][0]+=0;

  // 10
  //   A_[1][1]+=
  addElementNorm( 10, 
	      (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/3., weightNrm_ );

  // 11
  //   A_[1][2]+=0;

  // 12
  //   A_[1][3]+=t;
  //   A_[3][1]+=t;
  addElementNorm( 12, 
	      (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12., weightNrm_ );

  // 13
  //   A_[1][4]+=t;
  //   A_[4][1]+=t;
  addElementNorm( 13, 
	      (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12., weightNrm_ );

  // 14
  //   A_[1][5]+=0;

  // 15
  //   A_[1][6]+=0;

  // 16
  //   A_[1][7]+=t;
  //   A_[7][1]+=t;
  addElementNorm( 16, 
	      (y0 + y1 + y2)/3., weightNrm_ );

  // 17
  //   A_[1][8]+=0;

  // 18
  //   A_[1][9]+=0;

  //   A_[2][0]+=0;

  //   A_[2][1]+=0;

  // 19
  //   A_[2][2]+=
  addElementNorm( 19, 
	      (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/3., weightNrm_ );

  // 20
  //   A_[2][3]+=0;

  // 21
  //   A_[2][4]+=t;
  //   A_[4][2]+=t;
  addElementNorm( 21, 
	      (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12., weightNrm_ );

  // 22
  //   A_[2][5]+=t;
  //   A_[5][2]+=t;
  addElementNorm( 22, 
	      (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12., weightNrm_ );

  // 23
  //   A_[2][6]+=0;

  // 24
  //   A_[2][7]+=0;

  // 25
  //   A_[2][8]+=t;
  //   A_[8][2]+=t;
  addElementNorm( 25, 
	      (z0 + z1 + z2)/3., weightNrm_ );

  // 26
  //   A_[2][9]+=0;

  //   A_[3][0]+=(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12.;
  //   A_[3][1]+=(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12.;
  //   A_[3][2]+=0;

  // 27
  //   A_[3][3]+=
  addElementNorm( 27, 
	      (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12. + 
	      (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12., weightNrm_ );

  // 28
  //   A_[3][4]+=t;
  //   A_[4][3]+=t;
  addElementNorm( 28, 
	      (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24., weightNrm_ );

  // 29
  //   A_[3][5]+=t;
  //   A_[5][3]+=t;
  addElementNorm( 29, 
	      (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24., weightNrm_ );

  // 30
  //   A_[3][6]+=t;
  //   A_[6][3]+=t;
  addElementNorm( 30, 
	      (y0 + y1 + y2)/6., weightNrm_ );

  // 31
  //   A_[3][7]+=t;
  //   A_[7][3]+=t;
  addElementNorm( 31, 
	      (x0 + x1 + x2)/6., weightNrm_ );

  // 32
  //   A_[3][8]+=0;

  // 33
  //   A_[3][9]+=0;

  //   A_[4][0]+=0;
  //   A_[4][1]+=(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12.;
  //   A_[4][2]+=(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12.;
  //   A_[4][3]+=(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24.;

  // 34
  //   A_[4][4]+=
  addElementNorm( 34,
	      (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12. + 
	      (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12., weightNrm_ );

  // 35
  //   A_[4][5]+=t;
  //   A_[5][4]+=t;
  addElementNorm( 35, 
	      (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24., weightNrm_ );

  // 36
  //   A_[4][6]+=0;

  // 37
  //   A_[4][7]+=t;
  //   A_[7][4]+=t;
  addElementNorm( 37, 
	      (z0 + z1 + z2)/6., weightNrm_ );

  // 38
  //   A_[4][8]+=t;
  //   A_[8][4]+=t;
  addElementNorm( 38, 
	      (y0 + y1 + y2)/6., weightNrm_ );

  // 39
  //   A_[4][9]+=0;

  //   A_[5][0]+=(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12.;
  //   A_[5][1]+=0;
  //   A_[5][2]+=(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12.;
  //   A_[5][3]+=(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24.;
  //   A_[5][4]+=(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24.;

  // 40
  //   A_[5][5]+=
  addElementNorm( 40, 
	      (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12. + 
	      (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12., weightNrm_ );

  // 41
  //   A_[5][6]+=t;
  //   A_[6][5]+=t;
  addElementNorm( 41, 
	      (z0 + z1 + z2)/6., weightNrm_ );

  // 42
  //   A_[5][7]+=0;

  // 43
  //   A_[5][8]+=t;
  //   A_[8][5]+=t;
  addElementNorm( 43, 
	      (x0 + x1 + x2)/6., weightNrm_ );

  // 44
  //   A_[5][9]+=0;

  //   A_[6][0]+=(x0 + x1 + x2)/3.;
  //   A_[6][1]+=0;
  //   A_[6][2]+=0;
  //   A_[6][3]+=(y0 + y1 + y2)/6.;
  //   A_[6][4]+=0;
  //   A_[6][5]+=(z0 + z1 + z2)/6.;

  // 45
  //   A_[6][6]+=
  addElementNorm( 45, 
	      0.5, weightNrm_ );

  // 46
  //   A_[6][7]+=0;

  // 47
  //   A_[6][8]+=0;

  // 48
  //   A_[6][9]+=0;

  //   A_[7][0]+=0;
  //   A_[7][1]+=(y0 + y1 + y2)/3.;
  //   A_[7][2]+=0;
  //   A_[7][3]+=(x0 + x1 + x2)/6.;
  //   A_[7][4]+=(z0 + z1 + z2)/6.;
  //   A_[7][5]+=0;
  //   A_[7][6]+=0;

  // 49
  //   A_[7][7]+=
  addElementNorm( 49, 
	      0.5, weightNrm_ );

  // 50
  //   A_[7][8]+=0;

  // 51
  //   A_[7][9]+=0;

  //   A_[8][0]+=0;
  //   A_[8][1]+=0;
  //   A_[8][2]+=(z0 + z1 + z2)/3.;
  //   A_[8][3]+=0;
  //   A_[8][4]+=(y0 + y1 + y2)/6.;
  //   A_[8][5]+=(x0 + x1 + x2)/6.;
  //   A_[8][6]+=0;
  //   A_[8][7]+=0;

  // 52
  //   A_[8][8]+=
  addElementNorm( 52, 
	      0.5, weightNrm_ );

  // 53
  //   A_[8][9]+=0;

  //   A_[9][0]+=0;
  //   A_[9][1]+=0;
  //   A_[9][2]+=0;
  //   A_[9][3]+=0;
  //   A_[9][4]+=0;
  //   A_[9][5]+=0;
  //   A_[9][6]+=0;
  //   A_[9][7]+=0;
  //   A_[9][8]+=0;

  // 54
  //   A_[9][9]+=0;
}

void PIEM::addPIEMElementsNormB( Point3d& p0, Point3d& p1, Point3d& p2, Vector3d& nm )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double z0 = p0.z;
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;
  double x2 = p2.x;
  double y2 = p2.y;
  double z2 = p2.z;
  double nx = nm.x;
  double ny = nm.y;
  double nz = nm.z;

  // 55
//   b_[0] = 
  addElementNorm( 55, 
	      (nx*(x0 + x1 + x2))/3., weightNrm_ );

  // 56
//   b_[1] = 
  addElementNorm( 56, 
	      (ny*(y0 + y1 + y2))/3., weightNrm_ );

  // 57
//   b_[2] = 
  addElementNorm( 57, 
	      (nz*(z0 + z1 + z2))/3., weightNrm_ );

  // 58
//   b_[3] = 
  addElementNorm( 58, 
	      (ny*(x0 + x1 + x2))/6. + (nx*(y0 + y1 + y2))/6., weightNrm_ );

  // 59
//   b_[4] = 
  addElementNorm( 59, 
	      (nz*(y0 + y1 + y2))/6. + (ny*(z0 + z1 + z2))/6., weightNrm_ );

  // 60
//   b_[5] = 
  addElementNorm( 60, 
	      (nz*(x0 + x1 + x2))/6. + (nx*(z0 + z1 + z2))/6., weightNrm_ );

  // 61
//   b_[6] = 
  addElementNorm( 61, 
	      nx/2., weightNrm_ );

  // 62
//   b_[7] = 
  addElementNorm( 62, 
	      ny/2., weightNrm_ );

  // 63
//   b_[8] = 
  addElementNorm( 63, 
	      nz/2., weightNrm_ );

  // 64
//   b_[9] = 0;
}

void PIEM::addPIEMElementsNormC( Vector3d& nm )
{
  // 65
//   c_ = 
  addElementNorm( 65,
	      nm.x * nm.x + nm.y * nm.y + nm.z * nm.z, weightNrm_ );
}

double PIEM::error( std::vector<double>& coeff )
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

#if _USE_MKL_
bool PIEM::optimizeMKL( std::vector<double>& opt )
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

  NRSVD::svdcmp( A_, 10, 10, &w[0], v );

  // MKL SVD
  itpp::mat A(10,10);
//   mat A = randn(10,10);

  itpp::mat U, V;
  itpp::vec S;

#if 1
  for ( int i = 0; i < 10; ++i )
    {
      for ( int j = 0; j < 10; ++j )
	{
	  A.set(i, j, A_[i+1][j+1]);
	  cout << A_[i+1][j+1] << " ";
	}
      cout << endl;
      
    }
#endif

  cout << "A = " << A << endl;

//   cout << "U: " << endl;
//   cout << U << endl;
  cout << "S: " << endl;
  cout << S << endl;
//   cout << "V: " << endl;
//   cout << V << endl;

  cout << "***************************" << endl;
//   itpp::svd(A, S);
  itpp::svd(A, U, S, V);

//   cout << "U: " << endl;
//   cout << U << endl;
  cout << "S: " << endl;
  cout << S << endl;
//   cout << "V: " << endl;
//   cout << V << endl;

//   cout << "norm(A - U*diag(S)*V^T) = "
//       << itpp::round_to_zero(itpp::norm(A - U * itpp::diag(S) * itpp::transpose(V))) << endl;

  for ( int i = 0; i < 10; ++i )
    {
      cout << "i = " << i << " w = " << w[i+1] << " s = " << S[i] << endl;
    }

  double wmax=0.0f;
  for (int k = 1; k < NUM_VEC; ++k )
    {
      // 計算に失敗
#if defined(WIN32)
      //if ( std::_isnan( w[k] ) ) return false;
      if ( std::isnan( w[k] ) ) return false;
#else
      if ( isnan( w[k] ) ) return false;
#endif
      if ( std::fabs(w[k]) > wmax ) wmax = std::fabs(w[k]);
    }

  double wmin=wmax*0.000001f;
  for (int k = 1; k < NUM_VEC; ++k )
    if ( std::fabs(w[k]) < wmin ) w[k] = 0.;
  
  std::vector<double> x( NUM_VEC );
  NRSVD::svbksb( A_, &w[0], v, 10, 10, &b_[0], &x[0] );

  for ( int i = 0; i < NUM_VEC - 1; ++i )
    opt[i] = x[i+1];

  for ( int i = 1; i < NUM_VEC; ++i ) delete v[i];
  delete v;

  clearMatrixVec();

  return true;
}

#endif // _USE_MKL_

bool PIEM::optimize( std::vector<double>& opt )
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

  NRSVD::svdcmp( A_, 10, 10, &w[0], v );

  double wmax=0.0f;
  for (int k = 1; k < NUM_VEC; ++k )
    {
      // 計算に失敗
#if defined(WIN32)
      // if ( std::_isnan( w[k] ) ) return false;
      if ( std::isnan( w[k] ) ) return false;
#else
      if ( isnan( w[k] ) ) return false;
#endif
      if ( std::fabs(w[k]) > wmax ) wmax = std::fabs(w[k]);
    }

  double wmin=wmax*0.000001f;
  for (int k = 1; k < NUM_VEC; ++k )
    if ( std::fabs(w[k]) < wmin ) w[k] = 0.;
  
  std::vector<double> x( NUM_VEC );
  NRSVD::svbksb( A_, &w[0], v, 10, 10, &b_[0], &x[0] );

  for ( int i = 0; i < NUM_VEC - 1; ++i )
    opt[i] = x[i+1];

  for ( int i = 1; i < NUM_VEC; ++i ) delete v[i];
  delete v;

  clearMatrixVec();

  return true;
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

//
// PIEM から A, b, c を作成
//
void PIEM::copyPIEMToMatrixVec()
{
  int n = 0;
  double alpha = 1.0;
  for ( int i = 1; i < NUM_VEC; ++i )
    {
      for ( int j = i; j < NUM_VEC; ++j )
	{
	  double t = elmDis_[n] + alpha * area_ * elmNorm_[n];
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
      b_[i] = elmDis_[n] + alpha * area_ * elmNorm_[n];
      ++n;
    }
  
  c_ = elmDis_[n] + alpha * area_ * elmNorm_[n];
}

/////////////////////////////////////////////////////////////////////////////////////////////////////

#if 0

//
// set a Matrix of distance error metric
//
void PIEM::setDisA( Point3d& p0, Point3d& p1, Point3d& p2 )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double z0 = p0.z;
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;
  double x2 = p2.x;
  double y2 = p2.y;
  double z2 = p2.z;

  double p4x0 = x0*x0*x0*x0;
  double p3x0 = x0*x0*x0;
  double p2x0 = x0*x0;
  double p4x1 = x1*x1*x1*x1;
  double p3x1 = x1*x1*x1;
  double p2x1 = x1*x1;
  double p4x2 = x2*x2*x2*x2;
  double p3x2 = x2*x2*x2;
  double p2x2 = x2*x2;

  double p4y0 = y0*y0*y0*y0;
  double p3y0 = y0*y0*y0;
  double p2y0 = y0*y0;
  double p4y1 = y1*y1*y1*y1;
  double p3y1 = y1*y1*y1;
  double p2y1 = y1*y1;
  double p4y2 = y2*y2*y2*y2;
  double p3y2 = y2*y2*y2;
  double p2y2 = y2*y2;

  double p4z0 = z0*z0*z0*z0;
  double p3z0 = z0*z0*z0;
  double p2z0 = z0*z0;
  double p4z1 = z1*z1*z1*z1;
  double p3z1 = z1*z1*z1;
  double p2z1 = z1*z1;
  double p4z2 = z2*z2*z2*z2;
  double p3z2 = z2*z2*z2;
  double p2z2 = z2*z2;

  A_[0][0] = (p4x0 + p4x1 + p3x1*x2 + p2x1*p2x2 + 
	      x1*p3x2 + p4x2 + p3x0*(x1 + x2) + 
	      x0*(x1 + x2)*(p2x1 + p2x2) + 
	      p2x0*(p2x1 + x1*x2 + p2x2))/30.;

  A_[0][1] = A_[1][0] = (p2x2*(p2y0 + y0*y1 + p2y1 + 3*(y0 + y1)*y2 + 
			       6*p2y2) + x1*x2*
			 (p2y0 + 3*p2y1 + 4*y1*y2 + 3*p2y2 + 
			  2*y0*(y1 + y2)) + p2x0*
			 (6*p2y0 + p2y1 + y1*y2 + p2y2 + 3*y0*(y1 + y2)) + 
			 p2x1*(p2y0 + 6*p2y1 + 3*y1*y2 + p2y2 + 
			       y0*(3*y1 + y2)) 
			 + x0*(x1* (3*p2y0 + 4*y0*y1 + 3*p2y1 + 2*(y0 + y1)*y2 + 
				    p2y2) + x2*(3*p2y0 + p2y1 + 2*y1*y2 + 
						3*p2y2 + 2*y0*(y1 + 2*y2))))/180.;

  A_[0][2] = A_[2][0] = (p2x2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
			       6*p2z2) + x1*x2*
			 (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
			  2*z0*(z1 + z2)) + p2x0*
			 (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
			 p2x1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
			       z0*(3*z1 + z2)) 
			 + x0*(x1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
				    p2z2) + x2*(3*p2z0 + p2z1 + 2*z1*z2 + 
						3*p2z2 + 2*z0*(z1 + 2*z2))))/180.;

  A_[0][3] = A_[3][0] = (p3x0*(4*y0 + y1 + y2) + p3x1*(y0 + 4*y1 + y2) + 
			 p2x1*x2*(y0 + 3*y1 + 2*y2) + x1*p2x2*(y0 + 2*y1 + 3*y2) + 
			 p3x2*(y0 + y1 + 4*y2) + 
			 p2x0*(x1*(3*y0 + 2*y1 + y2) + x2*(3*y0 + y1 + 2*y2)) + 
			 x0*(2*x1*x2*(y0 + y1 + y2) + p2x1*(2*y0 + 3*y1 + y2) + 
			     p2x2*(2*y0 + y1 + 3*y2)))/120.;

  A_[0][4] = A_[4][0] = (p2x0*(3*y0*(4*z0 + z1 + z2) + y1*(3*z0 + 2*z1 + z2) + 
			       y2*(3*z0 + z1 + 2*z2)) + 
			 p2x1*(y0*(2*z0 + 3*z1 + z2) + 3*y1*(z0 + 4*z1 + z2) + 
			       y2*(z0 + 3*z1 + 2*z2)) + 
			 2*x1*x2*(y0*(z0 + z1 + z2) + y1*(z0 + 3*z1 + 2*z2) + 
				  y2*(z0 + 2*z1 + 3*z2)) + 
			 p2x2*(y0*(2*z0 + z1 + 3*z2) + y1*(z0 + 2*z1 + 3*z2) + 
			       3*y2*(z0 + z1 + 4*z2)) + 
			 2*x0*(x1*(y2*(z0 + z1 + z2) + y0*(3*z0 + 2*z1 + z2) + 
				   y1*(2*z0 + 3*z1 + z2)) + 
			       x2*(y1*(z0 + z1 + z2) + y0*(3*z0 + z1 + 2*z2) + 
				   y2*(2*z0 + z1 + 3*z2))))/360.;

  A_[0][5] = A_[5][0] = (p3x0*(4*z0 + z1 + z2) + p3x1*(z0 + 4*z1 + z2) + 
			 p2x1*x2*(z0 + 3*z1 + 2*z2) + x1*p2x2*(z0 + 2*z1 + 3*z2) + 
			 p3x2*(z0 + z1 + 4*z2) + 
			 p2x0*(x1*(3*z0 + 2*z1 + z2) + x2*(3*z0 + z1 + 2*z2)) + 
			 x0*(2*x1*x2*(z0 + z1 + z2) + p2x1*(2*z0 + 3*z1 + z2) + 
			     p2x2*(2*z0 + z1 + 3*z2)))/120.;

  A_[0][6] = A_[6][0] = (p3x0 + p2x0*(x1 + x2) + 
			 (x1 + x2)*(p2x1 + p2x2) + 
			 x0*(p2x1 + x1*x2 + p2x2))/20.;
  
  A_[0][7] = A_[7][0] = (p2x0*(3*y0 + y1 + y2) + p2x1*(y0 + 3*y1 + y2) + 
			 p2x2*(y0 + y1 + 3*y2) + x1*x2*(y0 + 2*(y1 + y2)) + 
			 x0*(x1*(2*y0 + 2*y1 + y2) + x2*(2*y0 + y1 + 2*y2)))/60.;

  A_[0][8] = A_[8][0] = (p2x0*(3*z0 + z1 + z2) + p2x1*(z0 + 3*z1 + z2) + 
			 p2x2*(z0 + z1 + 3*z2) + x1*x2*(z0 + 2*(z1 + z2)) + 
			 x0*(x1*(2*z0 + 2*z1 + z2) + x2*(2*z0 + z1 + 2*z2)))/60.;

  A_[0][9] = A_[9][0] = (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12.;


  A_[1][1] = (p4y0 + p4y1 + p3y1*y2 + p2y1*p2y2 + 
	      y1*p3y2 + p4y2 + p3y0*(y1 + y2) + 
	      y0*(y1 + y2)*(p2y1 + p2y2) + 
	      p2y0*(p2y1 + y1*y2 + p2y2))/30.;

  A_[1][2] = A_[2][1] = (p2y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
			       6*p2z2) + y1*y2*
			 (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
			  2*z0*(z1 + z2)) + p2y0*
			 (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
			 p2y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
			       z0*(3*z1 + z2)) 
			 + y0*(y1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
				    p2z2) + y2*(3*p2z0 + p2z1 + 2*z1*z2 + 
						3*p2z2 + 2*z0*(z1 + 2*z2))))/180.;

  A_[1][3] = A_[3][1] = (x2*(p3y0 + p2y0*y1 + y0*p2y1 + p3y1 + 
			     2*(p2y0 + y0*y1 + p2y1)*y2 + 
			     3*(y0 + y1)*p2y2 + 4*p3y2) + 
			 x0*(4*p3y0 + 3*p2y0*(y1 + y2) + 
			     (y1 + y2)*(p2y1 + p2y2) + 
			     2*y0*(p2y1 + y1*y2 + p2y2)) + 
			 x1*(p3y0 + 4*p3y1 + 3*p2y1*y2 + 2*y1*p2y2 + 
			     p3y2 + p2y0*(2*y1 + y2) + 
			     y0*(3*p2y1 + 2*y1*y2 + p2y2)))/120.;

  A_[1][4] = A_[4][1] = (p3y0*(4*z0 + z1 + z2) + p3y1*(z0 + 4*z1 + z2) + 
			 p2y1*y2*(z0 + 3*z1 + 2*z2) + y1*p2y2*(z0 + 2*z1 + 3*z2) + 
			 p3y2*(z0 + z1 + 4*z2) + 
			 p2y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2)) + 
			 y0*(2*y1*y2*(z0 + z1 + z2) + p2y1*(2*z0 + 3*z1 + z2) + 
			     p2y2*(2*z0 + z1 + 3*z2)))/120.;

  A_[1][5] = A_[5][1] = (x1*(p2y0*(3*z0 + 2*z1 + z2) + 3*p2y1*(z0 + 4*z1 + z2) + 
			     2*y1*y2*(z0 + 3*z1 + 2*z2) + p2y2*(z0 + 2*z1 + 3*z2) + 
			     2*y0*(y2*(z0 + z1 + z2) + y1*(2*z0 + 3*z1 + z2))) + 
			 x0*(2*y1*y2*(z0 + z1 + z2) + 3*p2y0*(4*z0 + z1 + z2) + 
			     p2y1*(2*z0 + 3*z1 + z2) + p2y2*(2*z0 + z1 + 3*z2) + 
			     2*y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2))) + 
			 x2*(p2y0*(3*z0 + z1 + 2*z2) + p2y1*(z0 + 3*z1 + 2*z2) + 
			     2*y1*y2*(z0 + 2*z1 + 3*z2) + 3*p2y2*(z0 + z1 + 4*z2) + 
			     2*y0*(y1*(z0 + z1 + z2) + y2*(2*z0 + z1 + 3*z2))))/360.;

  A_[1][6] = A_[6][1] = (x1*(p2y0 + 2*y0*y1 + 3*p2y1 + y0*y2 + 2*y1*y2 + 
			     p2y2) + x2*(p2y0 + y0*y1 + p2y1 + 2*y0*y2 + 
					 2*y1*y2 + 3*p2y2) + 
			 x0*(3*p2y0 + p2y1 + y1*y2 + p2y2 + 2*y0*(y1 + y2)))/60.;

  A_[1][7] = A_[7][1] = (p3y0 + p2y0*(y1 + y2) + 
			 (y1 + y2)*(p2y1 + p2y2) + 
			 y0*(p2y1 + y1*y2 + p2y2))/20.;

  A_[1][8] = A_[8][1] = (p2y0*(3*z0 + z1 + z2) + p2y1*(z0 + 3*z1 + z2) + 
			 p2y2*(z0 + z1 + 3*z2) + y1*y2*(z0 + 2*(z1 + z2)) + 
			 y0*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2)))/60.;

  A_[1][9] = A_[9][1] = (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12.;


  A_[2][2] = (p4z0 + p4z1 + p3z1*z2 + p2z1*p2z2 + 
	      z1*p3z2 + p4z2 + p3z0*(z1 + z2) + 
	      z0*(z1 + z2)*(p2z1 + p2z2) + 
	      p2z0*(p2z1 + z1*z2 + p2z2))/30.;

  A_[2][3] = A_[3][2] = (x1*(y0*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
				 p2z2) + y2*(p2z0 + 3*p2z1 + 4*z1*z2 + 
					     3*p2z2 + 2*z0*(z1 + z2)) + 
			     2*y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
				   z0*(3*z1 + z2))) + x2*
			 (2*y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
				6*p2z2) + y1*
			  (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
			   2*z0*(z1 + z2)) + y0*
			  (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
			   2*z0*(z1 + 2*z2))) + 
			 x0*(y1*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
				 p2z2) + 2*y0*
			     (6*p2z0 + p2z1 + z1*z2 + p2z2 + 
			      3*z0*(z1 + z2)) + y2*
			     (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
			      2*z0*(z1 + 2*z2))))/360.;

  A_[2][4] = A_[4][2] = (y2*(p3z0 + p2z0*z1 + z0*p2z1 + p3z1 + 
			     2*(p2z0 + z0*z1 + p2z1)*z2 + 
			     3*(z0 + z1)*p2z2 + 4*p3z2) + 
			 y0*(4*p3z0 + 3*p2z0*(z1 + z2) + 
			     (z1 + z2)*(p2z1 + p2z2) + 
			     2*z0*(p2z1 + z1*z2 + p2z2)) + 
			 y1*(p3z0 + 4*p3z1 + 3*p2z1*z2 + 2*z1*p2z2 + 
			     p3z2 + p2z0*(2*z1 + z2) + 
			     z0*(3*p2z1 + 2*z1*z2 + p2z2)))/120.;

  A_[2][5] = A_[5][2] = (x2*(p3z0 + p2z0*z1 + z0*p2z1 + p3z1 + 
			     2*(p2z0 + z0*z1 + p2z1)*z2 + 
			     3*(z0 + z1)*p2z2 + 4*p3z2) + 
			 x0*(4*p3z0 + 3*p2z0*(z1 + z2) + 
			     (z1 + z2)*(p2z1 + p2z2) + 
			     2*z0*(p2z1 + z1*z2 + p2z2)) + 
			 x1*(p3z0 + 4*p3z1 + 3*p2z1*z2 + 2*z1*p2z2 + 
			     p3z2 + p2z0*(2*z1 + z2) + 
			     z0*(3*p2z1 + 2*z1*z2 + p2z2)))/120.;

  A_[2][6] = A_[6][2] = (x1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
			     p2z2) + x2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
					 2*z1*z2 + 3*p2z2) + 
			 x0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[2][7] = A_[7][2] = (y1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
			     p2z2) + y2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
					 2*z1*z2 + 3*p2z2) + 
			 y0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[2][8] = A_[8][2] = (p3z0 + p2z0*(z1 + z2) + 
			 (z1 + z2)*(p2z1 + p2z2) + 
			 z0*(p2z1 + z1*z2 + p2z2))/20.;

  A_[2][9] = A_[9][2] = (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12.;


  A_[3][3] = (p2x2*(p2y0 + y0*y1 + p2y1 + 3*(y0 + y1)*y2 + 
		    6*p2y2) + x1*x2*
	      (p2y0 + 3*p2y1 + 4*y1*y2 + 3*p2y2 + 
	       2*y0*(y1 + y2)) + p2x0*
	      (6*p2y0 + p2y1 + y1*y2 + p2y2 + 3*y0*(y1 + y2)) + 
	      p2x1*(p2y0 + 6*p2y1 + 3*y1*y2 + p2y2 + 
		    y0*(3*y1 + y2)) 
	      + x0*(x1* (3*p2y0 + 4*y0*y1 + 3*p2y1 + 2*(y0 + y1)*y2 + 
			 p2y2) + x2*(3*p2y0 + p2y1 + 2*y1*y2 + 
				     3*p2y2 + 2*y0*(y1 + 2*y2))))/180.;

  A_[3][4] = A_[4][3] = (x1*(p2y0*(3*z0 + 2*z1 + z2) + 3*p2y1*(z0 + 4*z1 + z2) + 
			     2*y1*y2*(z0 + 3*z1 + 2*z2) + p2y2*(z0 + 2*z1 + 3*z2) + 
			     2*y0*(y2*(z0 + z1 + z2) + y1*(2*z0 + 3*z1 + z2))) + 
			 x0*(2*y1*y2*(z0 + z1 + z2) + 3*p2y0*(4*z0 + z1 + z2) + 
			     p2y1*(2*z0 + 3*z1 + z2) + p2y2*(2*z0 + z1 + 3*z2) + 
			     2*y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2))) + 
			 x2*(p2y0*(3*z0 + z1 + 2*z2) + p2y1*(z0 + 3*z1 + 2*z2) + 
			     2*y1*y2*(z0 + 2*z1 + 3*z2) + 3*p2y2*(z0 + z1 + 4*z2) + 
			     2*y0*(y1*(z0 + z1 + z2) + y2*(2*z0 + z1 + 3*z2))))/360.;

  A_[3][5] = A_[5][3] = (p2x0*(3*y0*(4*z0 + z1 + z2) + y1*(3*z0 + 2*z1 + z2) + 
			       y2*(3*z0 + z1 + 2*z2)) + 
			 p2x1*(y0*(2*z0 + 3*z1 + z2) + 3*y1*(z0 + 4*z1 + z2) + 
			       y2*(z0 + 3*z1 + 2*z2)) + 
			 2*x1*x2*(y0*(z0 + z1 + z2) + y1*(z0 + 3*z1 + 2*z2) + 
				  y2*(z0 + 2*z1 + 3*z2)) + 
			 p2x2*(y0*(2*z0 + z1 + 3*z2) + y1*(z0 + 2*z1 + 3*z2) + 
			       3*y2*(z0 + z1 + 4*z2)) + 
			 2*x0*(x1*(y2*(z0 + z1 + z2) + y0*(3*z0 + 2*z1 + z2) + 
				   y1*(2*z0 + 3*z1 + z2)) + 
			       x2*(y1*(z0 + z1 + z2) + y0*(3*z0 + z1 + 2*z2) + 
				   y2*(2*z0 + z1 + 3*z2))))/360.;

  A_[3][6] = A_[6][3] = (p2x0*(3*y0 + y1 + y2) + p2x1*(y0 + 3*y1 + y2) + 
			 p2x2*(y0 + y1 + 3*y2) + x1*x2*(y0 + 2*(y1 + y2)) + 
			 x0*(x1*(2*y0 + 2*y1 + y2) + x2*(2*y0 + y1 + 2*y2)))/60.;

  A_[3][7] = A_[7][3] = (x1*(p2y0 + 2*y0*y1 + 3*p2y1 + y0*y2 + 2*y1*y2 + 
			     p2y2) + x2*(p2y0 + y0*y1 + p2y1 + 2*y0*y2 + 
					 2*y1*y2 + 3*p2y2) + 
			 x0*(3*p2y0 + p2y1 + y1*y2 + p2y2 + 2*y0*(y1 + y2)))/60.;

  A_[3][8] = A_[8][3] = (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
			     y2*(2*z0 + z1 + 2*z2)) + 
			 x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
			     y1*(z0 + 2*(z1 + z2))) + 
			 x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
			     y2*(z0 + 2*(z1 + z2))))/120.;

  A_[3][9] = A_[9][3] = (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24.;


  A_[4][4] = (p2y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		    6*p2z2) + y1*y2*
	      (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
	       2*z0*(z1 + z2)) + p2y0*
	      (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	      p2y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		    z0*(3*z1 + z2)) 
	      + y0*(y1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			 p2z2) + y2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				     3*p2z2 + 2*z0*(z1 + 2*z2))))/180.;

  A_[4][5] = A_[5][4] = (x1*(y0*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
				 p2z2) + y2*(p2z0 + 3*p2z1 + 4*z1*z2 + 
					     3*p2z2 + 2*z0*(z1 + z2)) + 
			     2*y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
				   z0*(3*z1 + z2))) + x2*
			 (2*y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
				6*p2z2) + y1*
			  (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
			   2*z0*(z1 + z2)) + y0*
			  (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
			   2*z0*(z1 + 2*z2))) + 
			 x0*(y1*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
				 p2z2) + 2*y0*
			     (6*p2z0 + p2z1 + z1*z2 + p2z2 + 
			      3*z0*(z1 + z2)) + y2*
			     (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
			      2*z0*(z1 + 2*z2))))/360.;

  A_[4][6] = A_[6][4] = (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
			     y2*(2*z0 + z1 + 2*z2)) + 
			 x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
			     y1*(z0 + 2*(z1 + z2))) + 
			 x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
			     y2*(z0 + 2*(z1 + z2))))/120.;

  A_[4][7] = A_[7][4] = (p2y0*(3*z0 + z1 + z2) + p2y1*(z0 + 3*z1 + z2) + 
			 p2y2*(z0 + z1 + 3*z2) + y1*y2*(z0 + 2*(z1 + z2)) + 
			 y0*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2)))/60.;

  A_[4][8] = A_[8][4] = (y1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
			     p2z2) + y2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
					 2*z1*z2 + 3*p2z2) + 
			 y0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[4][9] = A_[9][4] = (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24.;


  A_[5][5] = (p2x2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		    6*p2z2) + x1*x2*
	      (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
	       2*z0*(z1 + z2)) + p2x0*
	      (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	      p2x1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		    z0*(3*z1 + z2)) 
	      + x0*(x1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			 p2z2) + x2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				     3*p2z2 + 2*z0*(z1 + 2*z2))))/180.;

  A_[5][6] = A_[6][5] = (p2x0*(3*z0 + z1 + z2) + p2x1*(z0 + 3*z1 + z2) + 
			 p2x2*(z0 + z1 + 3*z2) + x1*x2*(z0 + 2*(z1 + z2)) + 
			 x0*(x1*(2*z0 + 2*z1 + z2) + x2*(2*z0 + z1 + 2*z2)))/60.;

  A_[5][7] = A_[7][5] = (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
			     y2*(2*z0 + z1 + 2*z2)) + 
			 x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
			     y1*(z0 + 2*(z1 + z2))) + 
			 x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
			     y2*(z0 + 2*(z1 + z2))))/120.;

  A_[5][8] = A_[8][5] = (x1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
			     p2z2) + x2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
					 2*z1*z2 + 3*p2z2) + 
			 x0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[5][9] = A_[9][5] = (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24.;


  A_[6][6] = (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12.;

  A_[6][7] = A_[7][6] = (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24.;

  A_[6][8] = A_[8][6] = (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24.;

  A_[6][9] = A_[9][6] = (x0 + x1 + x2)/6.;


  A_[7][7] = (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12.;

  A_[7][8] = A_[8][7] = (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24.;

  A_[7][9] = A_[9][7] = (y0 + y1 + y2)/6.;


  A_[8][8] = (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12.;

  A_[8][9] = A_[9][8] = (z0 + z1 + z2)/6.;


  A_[9][9] = 0.5;

  ////////////////////////////////////////////////////////////////////////////////////////////

#if 0
  A_[1][0] = (p2x2*(p2y0 + y0*y1 + p2y1 + 3*(y0 + y1)*y2 + 
		    6*p2y2) + x1*x2*
	      (p2y0 + 3*p2y1 + 4*y1*y2 + 3*p2y2 + 
	       2*y0*(y1 + y2)) + p2x0*
	      (6*p2y0 + p2y1 + y1*y2 + p2y2 + 3*y0*(y1 + y2)) + 
	      p2x1*(p2y0 + 6*p2y1 + 3*y1*y2 + p2y2 + 
		    y0*(3*y1 + y2)) 
	      + x0*(x1* (3*p2y0 + 4*y0*y1 + 3*p2y1 + 2*(y0 + y1)*y2 + 
			 p2y2) + x2*(3*p2y0 + p2y1 + 2*y1*y2 + 
				     3*p2y2 + 2*y0*(y1 + 2*y2))))/180.;
#endif

#if 0
  A_[2][0] = (p2x2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		    6*p2z2) + x1*x2*
	      (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
	       2*z0*(z1 + z2)) + p2x0*
	      (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	      p2x1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		    z0*(3*z1 + z2)) 
	      + x0*(x1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			 p2z2) + x2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				     3*p2z2 + 2*z0*(z1 + 2*z2))))/180.;

  A_[2][1] = (p2y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		    6*p2z2) + y1*y2*
	      (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
	       2*z0*(z1 + z2)) + p2y0*
	      (6*p2z0 + p2z1 + z1*z2 + p2z2 + 3*z0*(z1 + z2)) + 
	      p2y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
		    z0*(3*z1 + z2)) 
	      + y0*(y1* (3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
			 p2z2) + y2*(3*p2z0 + p2z1 + 2*z1*z2 + 
				     3*p2z2 + 2*z0*(z1 + 2*z2))))/180.;
#endif

#if 0
  A_[3][0] = (p3x0*(4*y0 + y1 + y2) + p3x1*(y0 + 4*y1 + y2) + 
	      p2x1*x2*(y0 + 3*y1 + 2*y2) + x1*p2x2*(y0 + 2*y1 + 3*y2) + 
	      p3x2*(y0 + y1 + 4*y2) + 
	      p2x0*(x1*(3*y0 + 2*y1 + y2) + x2*(3*y0 + y1 + 2*y2)) + 
	      x0*(2*x1*x2*(y0 + y1 + y2) + p2x1*(2*y0 + 3*y1 + y2) + 
		  p2x2*(2*y0 + y1 + 3*y2)))/120.;

  A_[3][1] = (x2*(p3y0 + p2y0*y1 + y0*p2y1 + p3y1 + 
		  2*(p2y0 + y0*y1 + p2y1)*y2 + 
		  3*(y0 + y1)*p2y2 + 4*p3y2) + 
	      x0*(4*p3y0 + 3*p2y0*(y1 + y2) + 
		  (y1 + y2)*(p2y1 + p2y2) + 
		  2*y0*(p2y1 + y1*y2 + p2y2)) + 
	      x1*(p3y0 + 4*p3y1 + 3*p2y1*y2 + 2*y1*p2y2 + 
		  p3y2 + p2y0*(2*y1 + y2) + 
		  y0*(3*p2y1 + 2*y1*y2 + p2y2)))/120.;

  A_[3][2] = (x1*(y0*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		      p2z2) + y2*(p2z0 + 3*p2z1 + 4*z1*z2 + 
				  3*p2z2 + 2*z0*(z1 + z2)) + 
		  2*y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
			z0*(3*z1 + z2))) + x2*
	      (2*y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		     6*p2z2) + y1*
	       (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + z2)) + y0*
	       (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + 2*z2))) + 
	      x0*(y1*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		      p2z2) + 2*y0*
		  (6*p2z0 + p2z1 + z1*z2 + p2z2 + 
		   3*z0*(z1 + z2)) + y2*
		  (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		   2*z0*(z1 + 2*z2))))/360.;
#endif

#if 0
  A_[4][0] = (p2x0*(3*y0*(4*z0 + z1 + z2) + y1*(3*z0 + 2*z1 + z2) + 
		    y2*(3*z0 + z1 + 2*z2)) + 
	      p2x1*(y0*(2*z0 + 3*z1 + z2) + 3*y1*(z0 + 4*z1 + z2) + 
		    y2*(z0 + 3*z1 + 2*z2)) + 
	      2*x1*x2*(y0*(z0 + z1 + z2) + y1*(z0 + 3*z1 + 2*z2) + 
		       y2*(z0 + 2*z1 + 3*z2)) + 
	      p2x2*(y0*(2*z0 + z1 + 3*z2) + y1*(z0 + 2*z1 + 3*z2) + 
		    3*y2*(z0 + z1 + 4*z2)) + 
	      2*x0*(x1*(y2*(z0 + z1 + z2) + y0*(3*z0 + 2*z1 + z2) + 
			y1*(2*z0 + 3*z1 + z2)) + 
		    x2*(y1*(z0 + z1 + z2) + y0*(3*z0 + z1 + 2*z2) + 
			y2*(2*z0 + z1 + 3*z2))))/360.;

  A_[4][1] = (p3y0*(4*z0 + z1 + z2) + p3y1*(z0 + 4*z1 + z2) + 
	      p2y1*y2*(z0 + 3*z1 + 2*z2) + y1*p2y2*(z0 + 2*z1 + 3*z2) + 
	      p3y2*(z0 + z1 + 4*z2) + 
	      p2y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2)) + 
	      y0*(2*y1*y2*(z0 + z1 + z2) + p2y1*(2*z0 + 3*z1 + z2) + 
		  p2y2*(2*z0 + z1 + 3*z2)))/120.;

  A_[4][2] = (y2*(p3z0 + p2z0*z1 + z0*p2z1 + p3z1 + 
		  2*(p2z0 + z0*z1 + p2z1)*z2 + 
		  3*(z0 + z1)*p2z2 + 4*p3z2) + 
	      y0*(4*p3z0 + 3*p2z0*(z1 + z2) + 
		  (z1 + z2)*(p2z1 + p2z2) + 
		  2*z0*(p2z1 + z1*z2 + p2z2)) + 
	      y1*(p3z0 + 4*p3z1 + 3*p2z1*z2 + 2*z1*p2z2 + 
		  p3z2 + p2z0*(2*z1 + z2) + 
		  z0*(3*p2z1 + 2*z1*z2 + p2z2)))/120.;

  A_[4][3] = (x1*(p2y0*(3*z0 + 2*z1 + z2) + 3*p2y1*(z0 + 4*z1 + z2) + 
		  2*y1*y2*(z0 + 3*z1 + 2*z2) + p2y2*(z0 + 2*z1 + 3*z2) + 
		  2*y0*(y2*(z0 + z1 + z2) + y1*(2*z0 + 3*z1 + z2))) + 
	      x0*(2*y1*y2*(z0 + z1 + z2) + 3*p2y0*(4*z0 + z1 + z2) + 
		  p2y1*(2*z0 + 3*z1 + z2) + p2y2*(2*z0 + z1 + 3*z2) + 
		  2*y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2))) + 
	      x2*(p2y0*(3*z0 + z1 + 2*z2) + p2y1*(z0 + 3*z1 + 2*z2) + 
		  2*y1*y2*(z0 + 2*z1 + 3*z2) + 3*p2y2*(z0 + z1 + 4*z2) + 
		  2*y0*(y1*(z0 + z1 + z2) + y2*(2*z0 + z1 + 3*z2))))/360.;
#endif

#if 0
  A_[5][0] = (p3x0*(4*z0 + z1 + z2) + p3x1*(z0 + 4*z1 + z2) + 
	      p2x1*x2*(z0 + 3*z1 + 2*z2) + x1*p2x2*(z0 + 2*z1 + 3*z2) + 
	      p3x2*(z0 + z1 + 4*z2) + 
	      p2x0*(x1*(3*z0 + 2*z1 + z2) + x2*(3*z0 + z1 + 2*z2)) + 
	      x0*(2*x1*x2*(z0 + z1 + z2) + p2x1*(2*z0 + 3*z1 + z2) + 
		  p2x2*(2*z0 + z1 + 3*z2)))/120.;

  A_[5][1] = (x1*(p2y0*(3*z0 + 2*z1 + z2) + 3*p2y1*(z0 + 4*z1 + z2) + 
		  2*y1*y2*(z0 + 3*z1 + 2*z2) + p2y2*(z0 + 2*z1 + 3*z2) + 
		  2*y0*(y2*(z0 + z1 + z2) + y1*(2*z0 + 3*z1 + z2))) + 
	      x0*(2*y1*y2*(z0 + z1 + z2) + 3*p2y0*(4*z0 + z1 + z2) + 
		  p2y1*(2*z0 + 3*z1 + z2) + p2y2*(2*z0 + z1 + 3*z2) + 
		  2*y0*(y1*(3*z0 + 2*z1 + z2) + y2*(3*z0 + z1 + 2*z2))) + 
	      x2*(p2y0*(3*z0 + z1 + 2*z2) + p2y1*(z0 + 3*z1 + 2*z2) + 
		  2*y1*y2*(z0 + 2*z1 + 3*z2) + 3*p2y2*(z0 + z1 + 4*z2) + 
		  2*y0*(y1*(z0 + z1 + z2) + y2*(2*z0 + z1 + 3*z2))))/360.;

  A_[5][2] = (x2*(p3z0 + p2z0*z1 + z0*p2z1 + p3z1 + 
		  2*(p2z0 + z0*z1 + p2z1)*z2 + 
		  3*(z0 + z1)*p2z2 + 4*p3z2) + 
	      x0*(4*p3z0 + 3*p2z0*(z1 + z2) + 
		  (z1 + z2)*(p2z1 + p2z2) + 
		  2*z0*(p2z1 + z1*z2 + p2z2)) + 
	      x1*(p3z0 + 4*p3z1 + 3*p2z1*z2 + 2*z1*p2z2 + 
		  p3z2 + p2z0*(2*z1 + z2) + 
		  z0*(3*p2z1 + 2*z1*z2 + p2z2)))/120.;

  A_[5][3] = (p2x0*(3*y0*(4*z0 + z1 + z2) + y1*(3*z0 + 2*z1 + z2) + 
		    y2*(3*z0 + z1 + 2*z2)) + 
	      p2x1*(y0*(2*z0 + 3*z1 + z2) + 3*y1*(z0 + 4*z1 + z2) + 
		    y2*(z0 + 3*z1 + 2*z2)) + 
	      2*x1*x2*(y0*(z0 + z1 + z2) + y1*(z0 + 3*z1 + 2*z2) + 
		       y2*(z0 + 2*z1 + 3*z2)) + 
	      p2x2*(y0*(2*z0 + z1 + 3*z2) + y1*(z0 + 2*z1 + 3*z2) + 
		    3*y2*(z0 + z1 + 4*z2)) + 
	      2*x0*(x1*(y2*(z0 + z1 + z2) + y0*(3*z0 + 2*z1 + z2) + 
			y1*(2*z0 + 3*z1 + z2)) + 
		    x2*(y1*(z0 + z1 + z2) + y0*(3*z0 + z1 + 2*z2) + 
			y2*(2*z0 + z1 + 3*z2))))/360.;

  A_[5][4] = (x1*(y0*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		      p2z2) + y2*(p2z0 + 3*p2z1 + 4*z1*z2 + 
				  3*p2z2 + 2*z0*(z1 + z2)) + 
		  2*y1*(p2z0 + 6*p2z1 + 3*z1*z2 + p2z2 + 
			z0*(3*z1 + z2))) + x2*
	      (2*y2*(p2z0 + z0*z1 + p2z1 + 3*(z0 + z1)*z2 + 
		     6*p2z2) + y1*
	       (p2z0 + 3*p2z1 + 4*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + z2)) + y0*
	       (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		2*z0*(z1 + 2*z2))) + 
	      x0*(y1*(3*p2z0 + 4*z0*z1 + 3*p2z1 + 2*(z0 + z1)*z2 + 
		      p2z2) + 2*y0*
		  (6*p2z0 + p2z1 + z1*z2 + p2z2 + 
		   3*z0*(z1 + z2)) + y2*
		  (3*p2z0 + p2z1 + 2*z1*z2 + 3*p2z2 + 
		   2*z0*(z1 + 2*z2))))/360.;
#endif

#if 0
  A_[6][0] = (p3x0 + p2x0*(x1 + x2) + 
	      (x1 + x2)*(p2x1 + p2x2) + 
	      x0*(p2x1 + x1*x2 + p2x2))/20.;

  A_[6][1] = (x1*(p2y0 + 2*y0*y1 + 3*p2y1 + y0*y2 + 2*y1*y2 + 
		  p2y2) + x2*(p2y0 + y0*y1 + p2y1 + 2*y0*y2 + 
			      2*y1*y2 + 3*p2y2) + 
	      x0*(3*p2y0 + p2y1 + y1*y2 + p2y2 + 2*y0*(y1 + y2)))/60.;

  A_[6][2] = (x1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		  p2z2) + x2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			      2*z1*z2 + 3*p2z2) + 
	      x0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[6][3] = (p2x0*(3*y0 + y1 + y2) + p2x1*(y0 + 3*y1 + y2) + 
	      p2x2*(y0 + y1 + 3*y2) + x1*x2*(y0 + 2*(y1 + y2)) + 
	      x0*(x1*(2*y0 + 2*y1 + y2) + x2*(2*y0 + y1 + 2*y2)))/60.;

  A_[6][4] = (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
		  y2*(2*z0 + z1 + 2*z2)) + 
	      x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
		  y1*(z0 + 2*(z1 + z2))) + 
	      x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
		  y2*(z0 + 2*(z1 + z2))))/120.;

  A_[6][5] = (p2x0*(3*z0 + z1 + z2) + p2x1*(z0 + 3*z1 + z2) + 
	      p2x2*(z0 + z1 + 3*z2) + x1*x2*(z0 + 2*(z1 + z2)) + 
	      x0*(x1*(2*z0 + 2*z1 + z2) + x2*(2*z0 + z1 + 2*z2)))/60.;
#endif

#if 0
  A_[7][0] = (p2x0*(3*y0 + y1 + y2) + p2x1*(y0 + 3*y1 + y2) + 
	      p2x2*(y0 + y1 + 3*y2) + x1*x2*(y0 + 2*(y1 + y2)) + 
	      x0*(x1*(2*y0 + 2*y1 + y2) + x2*(2*y0 + y1 + 2*y2)))/60.;

  A_[7][1] = (p3y0 + p2y0*(y1 + y2) + 
	      (y1 + y2)*(p2y1 + p2y2) + 
	      y0*(p2y1 + y1*y2 + p2y2))/20.;

  A_[7][2] = (y1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		  p2z2) + y2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			      2*z1*z2 + 3*p2z2) + 
	      y0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[7][3] = (x1*(p2y0 + 2*y0*y1 + 3*p2y1 + y0*y2 + 2*y1*y2 + 
		  p2y2) + x2*(p2y0 + y0*y1 + p2y1 + 2*y0*y2 + 
			      2*y1*y2 + 3*p2y2) + 
	      x0*(3*p2y0 + p2y1 + y1*y2 + p2y2 + 2*y0*(y1 + y2)))/60.;

  A_[7][4] = (p2y0*(3*z0 + z1 + z2) + p2y1*(z0 + 3*z1 + z2) + 
	      p2y2*(z0 + z1 + 3*z2) + y1*y2*(z0 + 2*(z1 + z2)) + 
	      y0*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2)))/60.;

  A_[7][5] = (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
		  y2*(2*z0 + z1 + 2*z2)) + 
	      x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
		  y1*(z0 + 2*(z1 + z2))) + 
	      x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
		  y2*(z0 + 2*(z1 + z2))))/120.;

  A_[7][6] = (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24.;
#endif

#if 0
  A_[8][0] = (p2x0*(3*z0 + z1 + z2) + p2x1*(z0 + 3*z1 + z2) + 
	      p2x2*(z0 + z1 + 3*z2) + x1*x2*(z0 + 2*(z1 + z2)) + 
	      x0*(x1*(2*z0 + 2*z1 + z2) + x2*(2*z0 + z1 + 2*z2)))/60.;

  A_[8][1] = (p2y0*(3*z0 + z1 + z2) + p2y1*(z0 + 3*z1 + z2) + 
	      p2y2*(z0 + z1 + 3*z2) + y1*y2*(z0 + 2*(z1 + z2)) + 
	      y0*(y1*(2*z0 + 2*z1 + z2) + y2*(2*z0 + z1 + 2*z2)))/60.;

  A_[8][2] = (p3z0 + p2z0*(z1 + z2) + 
	      (z1 + z2)*(p2z1 + p2z2) + 
	      z0*(p2z1 + z1*z2 + p2z2))/20.;

  A_[8][3] = (x0*(2*y0*(3*z0 + z1 + z2) + y1*(2*z0 + 2*z1 + z2) + 
		  y2*(2*z0 + z1 + 2*z2)) + 
	      x2*(y0*(2*z0 + z1 + 2*z2) + 2*y2*(z0 + z1 + 3*z2) + 
		  y1*(z0 + 2*(z1 + z2))) + 
	      x1*(y0*(2*z0 + 2*z1 + z2) + 2*y1*(z0 + 3*z1 + z2) + 
		  y2*(z0 + 2*(z1 + z2))))/120.;

  A_[8][4] = (y1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		  p2z2) + y2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			      2*z1*z2 + 3*p2z2) + 
	      y0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[8][5] = (x1*(p2z0 + 2*z0*z1 + 3*p2z1 + z0*z2 + 2*z1*z2 + 
		  p2z2) + x2*(p2z0 + z0*z1 + p2z1 + 2*z0*z2 + 
			      2*z1*z2 + 3*p2z2) + 
	      x0*(3*p2z0 + p2z1 + z1*z2 + p2z2 + 2*z0*(z1 + z2)))/60.;

  A_[8][6] = (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24.;

  A_[8][7] = (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24.;
#endif

#if 0
  A_[9][0] = (p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12.;

  A_[9][1] = (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12.;

  A_[9][2] = (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12.;

  A_[9][3] = (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24.;

  A_[9][4] = (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24.;

  A_[9][5] = (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24.;

  A_[9][6] = (x0 + x1 + x2)/6.;

  A_[9][7] = (y0 + y1 + y2)/6.;

  A_[9][8] = (z0 + z1 + z2)/6.;
#endif

}

//
// set a Matrix of normal error metric
//
void PIEM::setNormA( Point3d& p0, Point3d& p1, Point3d& p2 )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double z0 = p0.z;
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;
  double x2 = p2.x;
  double y2 = p2.y;
  double z2 = p2.z;
  
  double p2x0 = x0*x0;
  double p2x1 = x1*x1;
  double p2x2 = x2*x2;

  double p2y0 = y0*y0;
  double p2y1 = y1*y1;
  double p2y2 = y2*y2;

  double p2z0 = z0*z0;
  double p2z1 = z1*z1;
  double p2z2 = z2*z2;

  double t;

  A_[0][0]+=(p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/3.;
//   A_[0][1]+=0;
//   A_[0][2]+=0;
  t = (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12.;
  A_[0][3]+=t;
  A_[3][0]+=t;
//   A_[0][4]+=0;
  t = (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12.;
  A_[0][5]+=t;
  A_[5][0]+=t;
  t = (x0 + x1 + x2)/3.;
  A_[0][6]+=t;
  A_[6][0]+=t;
//   A_[0][7]+=0;
//   A_[0][8]+=0;
//   A_[0][9]+=0;
//   A_[1][0]+=0;
  A_[1][1]+=(p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/3.;
//   A_[1][2]+=0;
  t = (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12.;
  A_[1][3]+=t;
  A_[3][1]+=t;
  t = (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12.;
  A_[1][4]+=t;
  A_[4][1]+=t;
//   A_[1][5]+=0;
//   A_[1][6]+=0;
  t = (y0 + y1 + y2)/3.;
  A_[1][7]+=t;
  A_[7][1]+=t;
//   A_[1][8]+=0;
//   A_[1][9]+=0;
//   A_[2][0]+=0;
//   A_[2][1]+=0;
  A_[2][2]+=(p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/3.;
//   A_[2][3]+=0;
  t = (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12.;
  A_[2][4]+=t;
  A_[4][2]+=t;
  t = (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12.;
  A_[2][5]+=t;
  A_[5][2]+=t;
//   A_[2][6]+=0;
//   A_[2][7]+=0;
  t = (z0 + z1 + z2)/3.;
  A_[2][8]+=t;
  A_[8][2]+=t;
//   A_[2][9]+=0;
//   A_[3][0]+=(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12.;
//   A_[3][1]+=(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/12.;
//   A_[3][2]+=0;
  A_[3][3]+=(p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12. + 
    (p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12.;
  t = (x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24.;
  A_[3][4]+=t;
  A_[4][3]+=t;
  t = (y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24.;
  A_[3][5]+=t;
  A_[5][3]+=t;
  t = (y0 + y1 + y2)/6.;
  A_[3][6]+=t;
  A_[6][3]+=t;
  t = (x0 + x1 + x2)/6.;
  A_[3][7]+=t;
  A_[7][3]+=t;
//   A_[3][8]+=0;
//   A_[3][9]+=0;
//   A_[4][0]+=0;
//   A_[4][1]+=(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12.;
//   A_[4][2]+=(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/12.;
//   A_[4][3]+=(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/24.;
  A_[4][4]+=(p2y0 + p2y1 + y1*y2 + p2y2 + y0*(y1 + y2))/12. + 
    (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12.;
  t = (x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24.;
  A_[4][5]+=t;
  A_[5][4]+=t;
//   A_[4][6]+=0;
  t = (z0 + z1 + z2)/6.;
  A_[4][7]+=t;
  A_[7][4]+=t;
  t = (y0 + y1 + y2)/6.;
  A_[4][8]+=t;
  A_[8][4]+=t;
//   A_[4][9]+=0;
//   A_[5][0]+=(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12.;
//   A_[5][1]+=0;
//   A_[5][2]+=(x0*(2*z0 + z1 + z2) + x1*(z0 + 2*z1 + z2) + x2*(z0 + z1 + 2*z2))/12.;
//   A_[5][3]+=(y0*(2*z0 + z1 + z2) + y1*(z0 + 2*z1 + z2) + y2*(z0 + z1 + 2*z2))/24.;
//   A_[5][4]+=(x0*(2*y0 + y1 + y2) + x1*(y0 + 2*y1 + y2) + x2*(y0 + y1 + 2*y2))/24.;
  A_[5][5]+=(p2x0 + p2x1 + x1*x2 + p2x2 + x0*(x1 + x2))/12. + 
    (p2z0 + p2z1 + z1*z2 + p2z2 + z0*(z1 + z2))/12.;
  t = (z0 + z1 + z2)/6.;
  A_[5][6]+=t;
  A_[6][5]+=t;
//   A_[5][7]+=0;
  t = (x0 + x1 + x2)/6.;
  A_[5][8]+=t;
  A_[8][5]+=t;
//   A_[5][9]+=0;
//   A_[6][0]+=(x0 + x1 + x2)/3.;
//   A_[6][1]+=0;
//   A_[6][2]+=0;
//   A_[6][3]+=(y0 + y1 + y2)/6.;
//   A_[6][4]+=0;
//   A_[6][5]+=(z0 + z1 + z2)/6.;
  A_[6][6]+=0.5;
//   A_[6][7]+=0;
//   A_[6][8]+=0;
//   A_[6][9]+=0;
//   A_[7][0]+=0;
//   A_[7][1]+=(y0 + y1 + y2)/3.;
//   A_[7][2]+=0;
//   A_[7][3]+=(x0 + x1 + x2)/6.;
//   A_[7][4]+=(z0 + z1 + z2)/6.;
//   A_[7][5]+=0;
//   A_[7][6]+=0;
  A_[7][7]+=0.5;
//   A_[7][8]+=0;
//   A_[7][9]+=0;
//   A_[8][0]+=0;
//   A_[8][1]+=0;
//   A_[8][2]+=(z0 + z1 + z2)/3.;
//   A_[8][3]+=0;
//   A_[8][4]+=(y0 + y1 + y2)/6.;
//   A_[8][5]+=(x0 + x1 + x2)/6.;
//   A_[8][6]+=0;
//   A_[8][7]+=0;
  A_[8][8]+=0.5;
//   A_[8][9]+=0;
//   A_[9][0]+=0;
//   A_[9][1]+=0;
//   A_[9][2]+=0;
//   A_[9][3]+=0;
//   A_[9][4]+=0;
//   A_[9][5]+=0;
//   A_[9][6]+=0;
//   A_[9][7]+=0;
//   A_[9][8]+=0;
//   A_[9][9]+=0;
}

//
// set a vector of normal error metric
//
void PIEM::setNormB( Point3d& p0, Point3d& p1, Point3d& p2, Vector3d& nm )
{
  double x0 = p0.x;
  double y0 = p0.y;
  double z0 = p0.z;
  double x1 = p1.x;
  double y1 = p1.y;
  double z1 = p1.z;
  double x2 = p2.x;
  double y2 = p2.y;
  double z2 = p2.z;
  double nx = nm.x;
  double ny = nm.y;
  double nz = nm.z;

  b_[0] = (nx*(x0 + x1 + x2))/3.;
  b_[1] = (ny*(y0 + y1 + y2))/3.;
  b_[2] = (nz*(z0 + z1 + z2))/3.;
  b_[3] = (ny*(x0 + x1 + x2))/6. + (nx*(y0 + y1 + y2))/6.;
  b_[4] = (nz*(y0 + y1 + y2))/6. + (ny*(z0 + z1 + z2))/6.;
  b_[5] = (nz*(x0 + x1 + x2))/6. + (nx*(z0 + z1 + z2))/6.;
  b_[6] = nx/2.;
  b_[7] = ny/2.;
  b_[8] = nz/2.;
  b_[9] = 0;
}

//
// set a scalar of normal error metric
//
void PIEM::setNormC( Vector3d& nm )
{
  c_ = nm.x * nm.x + nm.y * nm.y + nm.z * nm.z;
}

#endif
