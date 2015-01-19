////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "MeshR.hxx"
#include "SMFRIO.hxx"
#include "SlimTree.hxx"
#include "SlimTreeIO.hxx"

int main( int argc, char* argv[] )
{
  if ( argc != 4 )
    {
      std::cerr << "Usage: " << argv[0] << " in.slim2t in.smf dis" << std::endl;
      exit( -1 );
    }

  SlimTreed slimtree;
  SlimTreedIO slimtreeIO( slimtree );
  slimtreeIO.inputFromFile( argv[1] );

#if 1

  double error = std::atof( argv[2] );
//   meshR.normalize();
//   meshR.scale( 10.0 );

  int count_ball = slimtree.countBalls( error );
  cout << "error: " << error << " ball count: " << count_ball << endl;

  std::ifstream ifs; ifs.open( "venus.pn" );

  double sum_dis = 0.0;
  double max_dis = 0.0;
  int num_points = 0;
  int out_points = 0;

  int i = 0;
  std::string cline;
  while ( getline(ifs, cline, '\n') )
    {
      std::istringstream isstr( cline );
      std::string str; 
      double x, y, z;
      isstr >> x; isstr >> y; isstr >> z;

      double nx, ny, nz;
      isstr >> nx; isstr >> ny; isstr >> nz;
      
      Point3d p( x, y, z );
#if 0
      double dis = slimtree.distance( p, error, DIST_AVR );
#endif

#if 1
      Vector3d n( nx, ny, nz );
      double dis = slimtree.grad_distance( p, n, error, DIST_AVR );
#endif

      if ( dis < 0 ) 
	{
	  ++( out_points );
	  continue;
	}

      sum_dis += dis;

      if ( i )
	{
	  if ( max_dis < dis ) max_dis = dis;
	}
      else max_dis = dis;

      ++i;
      ++( num_points );
      
    }

  ifs.close();

  sum_dis /= (double) num_points;
  
  cout << "avarage distance: " << sum_dis << endl;
  cout << "maximum distance: " << max_dis << endl;
  cout << out_points << " points is outside the balls. " << endl;
  
#endif

#if 0

  MeshRd meshR;
  SMFRdIO rio_input( meshR );
  rio_input.inputFromFile( argv[2] );

  double error = std::atof( argv[3] );
//   meshR.normalize();
//   meshR.scale( 10.0 );

  int count_ball = slimtree.countBalls( error );
  cout << "error: " << error << " ball count: " << count_ball << endl;

  double sum_dis = 0.0;
  double max_dis = 0.0;
  std::vector<double>& points = meshR.points();
  int num_points = meshR.numPoints();
  for ( int i = 0; i < meshR.numPoints(); ++i )
    {
      Point3d p( points[nXYZ * i], points[nXYZ * i + 1], points[nXYZ * i + 2] );
      double dis = slimtree.distance( p, error, DIST_AVR );
      if ( dis < 0 ) 
	{
	  --num_points;
	  continue;
	}

      sum_dis += dis;

      if ( i )
	{
	  if ( max_dis < dis ) max_dis = dis;
	}
      else max_dis = dis;
    }

//   cout << sum_dis << endl;
//   cout << num_points << endl;
  sum_dis /= (double) num_points;
  
  cout << "avarage distance: " << sum_dis << endl;
  cout << "maximum distance: " << max_dis << endl;
  cout << meshR.numPoints() - num_points << " points is outside the balls. " << endl;
  
#endif

  return 0;
}
