////////////////////////////////////////////////////////////////////
//
// $Id: mesh2slim.cxx 2021/06/04 16:33:15 kanai Exp $
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "MeshR.hxx"
#include "SMFRIO.hxx"
#include "SlimTree.hxx"
#include "SlimTreeIO.hxx"

#include "SlimSimplify.hxx"

#if 0
#include "mt19937ar.h"
Mt19937ar mt;
#endif

int main( int argc, char* argv[] )
{
  if ( argc != 3 )
    {
      std::cerr << "Usage: " << argv[0] << " in.obj out.slim2t " << std::endl;
      exit( -1 );
    }

  MeshRd meshR;
  SMFRdIO rio_input( meshR );
  rio_input.inputFromFile( argv[1] );
  
#if 0
  double w = 0.02;
  std::vector<double>& points = meshR.points();
  for ( int i = 0; i < meshR.numPoints(); ++i )
    {
      points[nXYZ * i] += ( w * (mt.genrand_real1() - 0.5) );
      points[nXYZ * i + 1] += ( w * (mt.genrand_real1() - 0.5) );
      points[nXYZ * i + 2] += ( w * (mt.genrand_real1() - 0.5) );
    }

  rio_input.outputToFile( "tmp.smf" );
#endif

#if 1
  meshR.normalize();
  meshR.scale( 10.0 );
//   rio_input.outputToFile( "tmp.smf" );
#endif

#if 1
  SlimTreed slimtree;
  SlimSimplify slimsimp( meshR, slimtree );
  slimsimp.setConstructType( SLIM_LEVEL );
  slimsimp.setIsStorePIEM( false );
//   slimsimp.setIsStorePIEM( true );
  slimsimp.apply();
  slimsimp.clear();
#endif

#if 1
  SlimTreedIO slimtreeIO( slimtree );
  slimtreeIO.outputToFile( argv[2] );
#endif

#if 0
  //
  // ランダム点を読み込んで誤差測定
  //
  for ( int k = 0; k < 3; ++k )
    {
      double error;
      if ( k == 0 ) error = 0.01;
      else if ( k == 1 ) error = 0.001;
      else error = 0.0001;

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

	  // "v"
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

      cout << "num. random points " << num_points << endl;

      sum_dis /= (double) num_points;
  
      cout << "avarage distance: " << sum_dis << endl;
      cout << "maximum distance: " << max_dis << endl;
      cout << out_points << " points is outside the balls. " << endl;
    }

#endif

#if 0
  //
  // ランダム点発生
  //
#define RAN_NUM 10

  mt.init();

  FILE* fp = fopen("tmp.pnt", "w");

  // 法線の生成
  meshR.createFaceNormals();

  for ( int i = 0; i < meshR.numFaces(); ++i )
    {
      Point3d p0, p1, p2;
      meshR.getFacePoints( i, p0, p1, p2 );

      Vector3d nrm;
      meshR.fnormal( i, nrm );

      for ( int j = 0; j < RAN_NUM; ++j )
	{
	  double a = mt.genrand_real1();
	  double b = mt.genrand_real1();
	  double c = mt.genrand_real1();
	  double d = a + b + c;
	  a /= d;	
	  b /= d;
	  c /= d;
	  Point3d p( a * p0.x + b * p1.x + c * p2.x,
		     a * p0.y + b * p1.y + c * p2.y,
		     a * p0.z + b * p1.z + c * p2.z );
 	  fprintf( fp, "%lf %lf %lf %lf %lf %lf\n", p.x, p.y, p.z, nrm.x, nrm.y, nrm.z );
	}

    }

  fclose( fp );

#endif

#if 0
  //
  // メッシュ上の頂点を使って誤差測定
  //
  double error = 0.0001;

  int count_ball = slimtree.countBalls( error );
  cout << "error: " << error << " ball count: " << count_ball << endl;

  double sum_dis = 0.0;
  double max_dis = 0.0;
  std::vector<double>& points = meshR.points();
  int num_points = meshR.numPoints();
  for ( int i = 0; i < meshR.numPoints(); ++i )
    {
      Point3d p( points[nXYZ * i], points[nXYZ * i + 1], points[nXYZ * i + 2] );
      double dis = slimtree.distance( p, error );
//       cout << "i " << i << " dis " << dis << endl;
      if ( dis < 0 ) 
	{
	  --num_points;
	  continue;
	}

//       cout << p << " " << dis << endl;
      sum_dis += dis;

      if ( i )
	{
	  if ( max_dis < dis ) max_dis = dis;
	}
      else max_dis = dis;
    }

  sum_dis /= (double) num_points;
  
  cout << "avarage distance: " << sum_dis << endl;
  cout << "maximum distance: " << max_dis << endl;
  cout << meshR.numPoints() - num_points << " points is outside the balls. " << endl;

#endif  

  return 0;
}
