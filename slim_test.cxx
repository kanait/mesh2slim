////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "envDep.h"

#include "Slim.hxx"
#include "SlimIO.hxx"

int main( int argc, char* argv[] )
{
  Slimf slim;
  SlimIOf slimIO( slim );

  slimIO.inputFromFile( argv[1] );

  SlimIOf wslimIO( slim );
  wslimIO.outputToFile("tmp.slim2");

  return 0;
}
