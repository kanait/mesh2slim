////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include <cstdio>
#include <cmath>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "Point2.h"
#include "Vector2.h"
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#include "PIEM2D.hxx"
#include "PIEM2DHandler.hxx"

#include "PointFit2D.hxx"
#include "PointFit2DHandler.hxx"

//////////////////////////////////////////////////////////////////////////////////////

void pointDistance( std::vector<Point2d>& pnt, 
		    std::vector<double>& coeff )
{
  double avrE = 0.0;
  double maxE = 0.0;

  for ( int i = 0; i < pnt.size(); ++i )
    {
      Point2d& p = pnt[i];
      double f0 = ( coeff[0] * p.x * p.x + coeff[1] * p.y * p.y +
		    coeff[2] * p.x * p.y + coeff[3] * p.x + 
		    coeff[4] * p.y + coeff[5] );
      double f = std::fabs( f0 );
      
      double dfx = 2.0 * coeff[0] * p.x + coeff[2] * p.y + coeff[3];
      double dfy = 2.0 * coeff[1] * p.y + coeff[2] * p.x + coeff[4];
      double df = std::sqrt( dfx * dfx + dfy * dfy );

      double dis = f / df;

      avrE += dis;

      if ( i ) 
	{
	  if ( maxE < dis ) maxE = dis;
	}
      else
	maxE = dis;
    }

  avrE /= (double) pnt.size();

  std::cout << "average: " << avrE << " max " << maxE << std::endl;
}

//////////////////////////////////////////////////////////////////////////////////////

void calcNormal( const Point2d& p0, const Point2d& p1, Vector2d& n ) 
{
  n.set( p0.y - p1.y, p1.x - p0.x );
  n.normalize();
}

double calcLength( const Point2d& p0, const Point2d& p1 ) 
{
  return p0.distance( p1 );
}

// regular sampling
void calcPointsNormalsLength_NU( std::vector<Point2d>& points,
			      std::vector<Point2d>& pnt,
			      std::vector<Vector2d>& nrm,
			      int num_point, 
			      double* sumLength )
{
  // 長さの合計
  *sumLength = 0.0;

  std::vector<double> polyLength( points.size() - 1 );
  for ( int i = 0; i < points.size()-1; ++i )
    {
      Point2d& p0 = points[i];
      Point2d& p1 = points[i+1];
      double length = calcLength( p0, p1 );
      polyLength[i] = length;
      *sumLength += ( length * length );
    }

  // min_length に対して num_point の数の点を発生
  // あとはなるべく等間隔になるようにする
  
  // calc point and normal
  for ( int i = 0; i < points.size()-1; ++i )
    {
      Point2d& p0 = points[i];
      Point2d& p1 = points[i+1];

      Vector2d n; calcNormal( p0, p1, n );
      
      for ( int j = 0; j < num_point; ++j )
	{
	  double param = (double) j / (double) num_point;
	  
	  Point2d p; p.interpolate( p0, p1, param );
	  pnt.push_back( p );
	  nrm.push_back( n );
// 	  std::cout << p << " " << n << std::endl;
	}
      
      
      // 最後の点
      if ( i == points.size()-2 )
	{
	  pnt.push_back( points[i+1] );
	  nrm.push_back( n );
	}
    }
  std::cout << "points " << pnt.size() << " normal " << nrm.size() << std::endl;
}


// spatially uniform sampling
void calcPointsNormalsLength( std::vector<Point2d>& points,
			      std::vector<Point2d>& pnt,
			      std::vector<Vector2d>& nrm,
			      int num_point, // 最も短いエッジに点を何点発生させるか
			      double* sumLength )
{
  // 長さの合計
  *sumLength = 0.0;

  // 最小長さ
  double min_length;

  std::vector<double> polyLength( points.size() - 1 );
  for ( int i = 0; i < points.size()-1; ++i )
    {
      Point2d& p0 = points[i];
      Point2d& p1 = points[i+1];
      double length = calcLength( p0, p1 );
      polyLength[i] = length;
      *sumLength += ( length * length );
      
      if ( i ) 
	{
	  if ( min_length > length ) min_length = length;
	}
      else
	min_length = length;
	  
    }

  // min_length に対して num_point の数の点を発生
  // あとはなるべく等間隔になるようにする
  
  // calc point and normal
  for ( int i = 0; i < points.size()-1; ++i )
    {
      Point2d& p0 = points[i];
      Point2d& p1 = points[i+1];

      int num_p = (int) ( (double) num_point * polyLength[i] / min_length );

      Vector2d n; calcNormal( p0, p1, n );
      
      for ( int j = 0; j < num_p; ++j )
	{
	  double param = (double) j / (double) num_p;
	  
	  Point2d p; p.interpolate( p0, p1, param );
	  pnt.push_back( p );
	  nrm.push_back( n );
// 	  std::cout << p << " " << n << std::endl;
	}
      
      
      // 最後の点
      if ( i == points.size()-2 )
	{
	  pnt.push_back( points[i+1] );
	  nrm.push_back( n );
	}
    }
  
  std::cout << "points " << pnt.size() << " normal " << nrm.size() << std::endl;
}

bool openFile( char* filename, std::vector<Point2d>& points )
{
  std::ifstream ifs; ifs.open( filename );
  if ( !ifs.is_open() )
    {
      std::cerr << "Cannot open " << filename << std::endl;
      return false;
    }

  std::string cline;
  getline(ifs, cline, '\n');
  std::istringstream isstr0( cline );
  int n;
  isstr0 >> n;
  points.resize( n );

  int i = 0;
  while ( getline(ifs, cline, '\n') )
    {
      std::istringstream isstr( cline );
      double x, y;
      isstr >> x; isstr >> y;
      points[i].set( x, y );
       
      ++i;
    }
    
  ifs.close();

  return true;
}

int main( int argc, char* argv[] )
{
  if ( argc != 3 )
    {
      std::cerr << "Usage: " << argv[0] << " file.txt num" << std::endl;
      exit( -1 );
    }

  std::vector<Point2d> points;
  openFile( argv[1], points );
 
  PIEM2DHandler piem_handle;

  PointFit2DHandler pf_handle;
  std::vector<Point2d> pnt;
  std::vector<Vector2d> nrm;
  double length;

  int num = atoi( argv[2] );

  int num_point = 200;
  if ( num == 0 )
    {
      piem_handle.setIsSetSumLength( false );
      piem_handle.setIsSetLengthWeight( false );
    }
  else if ( num == 1 )
    {
      piem_handle.setIsSetLengthWeight( true );
      piem_handle.setIsSetSumLength( false );
    }
  else if ( num == 2 )
    {
      piem_handle.setIsSetLengthWeight( true );
      piem_handle.setIsSetSumLength( true );
      calcPointsNormalsLength( points, pnt, nrm, num_point, &length );
    }
  else // num == 3
    {
      // spatially uniform sampling
      calcPointsNormalsLength( points, pnt, nrm, num_point, &length );
      // regular sampling
//       calcPointsNormalsLength_NU( points, pnt, nrm, num_point, &length );
    }
  

  std::vector<double> coeff( NUM_VEC - 1 );
  
  if ( num < 3 )
    {
      Point2d center;
      double radius;
      piem_handle.optimize( points, center, &radius, coeff );
    }
  else
    {
      pf_handle.optimize( pnt, nrm, length, coeff );
    }

  std::cout << coeff[0] << "x^2 + " << coeff[1] << "y^2 + " 
	    << coeff[2] << "x y + " << coeff[3] << "x + "
	    << coeff[4] << "y + " << coeff[5] << " = 0 " << std::endl;

  pointDistance( pnt, coeff );

  if ( num == 3 )
    {
      FILE* fp = fopen("tmp.txt", "w");

      fprintf( fp, "{\n");
      for ( int i = 0; i < pnt.size(); ++i )
	{
	  fprintf( fp, "{%g,%g},\n",pnt[i].x, pnt[i].y );
	}
      fprintf( fp, "}\n");
      
      fclose( fp );
    }


  return 0;
}
