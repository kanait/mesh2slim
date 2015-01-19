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

double func( double x, double y ) 
{
  double a = -0.5;
  double b = 0.0;
  double c = -0.5;
  double d = 0.0;
  double e = 0.0;
  double f = 0.5;
  return (a * x * x + b * x * y + c * y * y + d * x + e * y + f);
}

int main( int argc, char* argv[] )
{
  if ( argc != 2 )
    {
      std::cerr << "Usage: " << argv[0] << " out.smf(vobj) " << std::endl;
      exit( -1 );
    }

  MeshRd meshR;

  double size = 1.0;
  int div = 50;

  // create points
  for ( int i = 0; i <= div; ++i )
    {
      double x = -size + 2.0 * size * (double) i / (double) div;
      for ( int j = 0; j <= div; ++j )
	{
	  double y = -size + 2.0 * size * (double) j / (double) div;
	  double z = func( x, y );
	  meshR.addPoint( x, y, z );
	}
    }

  // create indices
  for ( int i = 0; i < div; ++i )
    {
      for ( int j = 0; j < div; ++j )
	{
	  unsigned int id0 = (div+1) * i + j;
	  unsigned int id1 = (div+1) * i + (j+1);
	  unsigned int id2 = (div+1) * (i+1) + j;
	  unsigned int id3 = (div+1) * (i+1) + (j+1);

	  meshR.addIndex( id0 );
	  meshR.addIndex( id2 );
	  meshR.addIndex( id1 );
  
	  meshR.addIndex( id1 );
	  meshR.addIndex( id2 );
	  meshR.addIndex( id3 );
	}
    }

  SMFRdIO rio_output( meshR );
  rio_output.outputToFile( argv[1] );

  return 0;
}
