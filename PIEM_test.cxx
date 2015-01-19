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

#include "BasisFunction.hxx"

#include "PIEM.hxx"
#include "PIEMHandler.hxx"
#include "SlimIO.hxx"

#define NUM_PARAM_QUAD 10

int main( int argc, char* argv[] )
{
  if ( argc != 3 )
    {
      std::cerr << "Usage: " << argv[0] << " in.smf(vobj) out.slim2 " << std::endl;
      exit( -1 );
    }

  MeshRd meshR;
  SMFRdIO rio_input( meshR );
  rio_input.inputFromFile( argv[1] );

  meshR.createFaceNormals();

  Slimd slim;
  slim.setDegree( 2 );
  slim.setOriented( true );
  slim.setLevel( 1 );
  slim.setBallArraySize( 2 ); // 2 * level
  std::vector< SlimBall<double> >& sballs = slim.sbarray( 0 );
  sballs.resize( 1 );

  PIEMHandler qiem_handle( meshR );
  Point3d center;
  double radius;
  qiem_handle.calcCenterRadius( center, &radius );
  Quadraticd* quad = new Quadraticd;
  qiem_handle.optimizeBasisFunction( center, quad->coeff() );
  
  sballs[0].setCenter( center );
  sballs[0].setSupport( radius );
  sballs[0].addBasisFunction( *quad );
  sballs[0].setIsLeaf( true );

//   Quadraticd* quad = new Quadraticd;
  std::vector<double>& coeff = quad->coeff();
  for ( int i = 0; i < QUADRATIC_COEFF; ++i )
    cout << "param " << i << " " << coeff[i] << endl;

  SlimIOd slimio( slim );
  slimio.outputToFile( argv[2] );
  return 0;
}
