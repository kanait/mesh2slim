////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIMIO_HXX
#define _SLIMIO_HXX 1

#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
using namespace std;

#include "SlimBall.hxx"
#include "BasisFunction.hxx"
#include "Slim.hxx"

template <typename T>
class SlimIO {

public:

  SlimIO( Slim<T>& slim ) { setSlim( slim ); };
  ~SlimIO(){};

  void setSlim( Slim<T>& slim ) { slim_ = &slim; };
  Slim<T>& slim() { return *slim_; };

  bool inputFromFile( const char* const infile ) {
#if 0
    ifstream in( infile, ios::in | ios::binary );
    if( !in ){ std::cerr << "cannot read " << infile << std::endl; return false; }
#endif
    
    FILE* in = std::fopen( infile, "r+b" );

    int N[2];
#if 0
    in.read( ( char* ) N, sizeof( N ) );
#endif
    std::fread( N, sizeof(int), 2, in );

    bool oriented;
    int degree = N[0];
    int level  = N[1];
    if( degree < 0 )
      {
	degree = -degree;
	oriented = false;
      }
    else
      oriented = true;
  
    slim().setDegree( degree );
    slim().setOriented( oriented );
    slim().setLevel( level );

    slim().setBallArraySize( 2 * level );

    for ( int i = 0; i < 2 * level; ++i )
      {
	std::vector< SlimBall<T> >& sballs = slim().sbarray(i);
	bool isLeaf = ( i%2 == 0 ) ? true : false;

	int n_ball;
#if 0
	in.read( (char*) &n_ball, sizeof(n_ball) );
#endif
	std::fread( &n_ball, sizeof(int), 1, in );
	if ( !n_ball ) continue;
	sballs.resize( n_ball );

	cout << "i = " << i << " num balls: " << n_ball << " isleaf " << isLeaf << endl;

	for ( int j = 0; j < n_ball; ++j )
	  {
	    sballs[j].setIsLeaf( isLeaf );

	    if ( degree == CUBIC )
	      {
		std::vector<T> c( CUBIC_COEFF_ALL );
#if 0
		in.read( (char*) &c[0], sizeof( c ) );
#endif
		std::fread( &c[0], sizeof(T), CUBIC_COEFF_ALL, in );

		sballs[j].setCenter( c[0], c[1], c[2] );
		sballs[j].setSupport( c[23] );
	      
		Cubic<T>* bf = new Cubic<T>;
		bf->setCoeffs( c );
		sballs[j].addBasisFunction( *bf );
	      }
	    else if ( degree == QUADRATIC )
	      {
		std::vector<T> c( QUADRATIC_COEFF_ALL );
#if 0
		in.read( (char*) &c[0], sizeof( c ) );
#endif
		std::fread(&c[0], sizeof(T), QUADRATIC_COEFF_ALL, in);

		sballs[j].setCenter( c[0], c[1], c[2] );
		sballs[j].setSupport( c[13] );

		Quadratic<T>* bf = new Quadratic<T>;
		bf->setCoeffs( c );
		sballs[j].addBasisFunction( *bf );

		cout << "\tballs: "; 
		for ( int k = 0; k < QUADRATIC_COEFF_ALL; ++k )
		  cout << c[k] << " ";
		cout << endl;
	      }
	    else
	      { 
		std::vector<T> c( LINEAR_COEFF_ALL );
#if 0
		in.read( (char*) &c[0], sizeof( c ) );
#endif
		std::fread( &c[0], sizeof(T), LINEAR_COEFF_ALL, in );

		sballs[j].setCenter( c[0], c[1], c[2] );
		sballs[j].setSupport( c[7] );

		Linear<T>* bf = new Linear<T>;
		bf->setCoeffs( c );
		sballs[j].addBasisFunction( *bf );
	      }
	  }
      }

#if 0  
    in.close();
#endif
    std::fclose( in );

    slim().print();

    return true;
  };

  bool outputToFile( const char* const outfile ) {
#if 1
    FILE* out = fopen( outfile, "w+b" );
#endif
#if 0
    ofstream out( outfile, ios::out | ios::binary );
    if( !out ) { std::cerr << "cannot read " << outfile << std::endl; return false; }
#endif
  
    int L[2];
    L[0] = ( slim().oriented() ) ? slim().degree() : -slim().degree();
    L[1] = slim().level();
#if 0
    out.write( (char*) L, sizeof( L ) );
#endif
#if 1
    fwrite( L, sizeof(int), 2, out );
#endif
    
    for( int i = 0; i < 2 * slim().level(); i++ )
      {
	std::vector< SlimBall<T> >& sballs = slim().sbarray(i);
#if 0
	int N = sballs.size();
	out.write( (char*) &N, sizeof( int ) );
#endif
#if 1
	L[0] = sballs.size();
	fwrite( L, sizeof(int), 1, out );
#endif

	for ( int j = 0; j < sballs.size(); ++j )
	  {
	    if ( slim().degree() == CUBIC )
	      {
		Cubic<T>& bf = (Cubic<T>&) sballs[j].bf();
		std::vector<T> c( CUBIC_COEFF_ALL );
		bf.getCoeffs( c );
		sballs[j].getCenter( &c[0] );
		sballs[j].getSupport( &(c[23]) );
#if 0
		out.write( (char*) &c[0], sizeof( c ) );
#endif
#if 1
		fwrite( &c[0], sizeof(float), CUBIC_COEFF_ALL, out );
#endif
	      }
	    else if ( slim().degree() == QUADRATIC )
	      {
		Quadratic<T>& bf = (Quadratic<T>&) sballs[j].bf();
		std::vector<T> c( QUADRATIC_COEFF_ALL );
		bf.getCoeffs( c );
		sballs[j].getCenter( &c[0] );
		sballs[j].getSupport( &(c[13]) );
		std::vector<float> cc( QUADRATIC_COEFF_ALL );
		for ( int i = 0; i < QUADRATIC_COEFF_ALL; ++i )
		  {
		    cc[i] = (float) c[i];
		    if ( std::fabs(cc[i]) < 1.0e-10 )
		      cc[i] = .0f;
		  }

#if 0
		out.write( (char*) &cc[0], sizeof( c ) );
#endif
#if 1
		fwrite( &cc[0], sizeof(float), QUADRATIC_COEFF_ALL, out );
#endif
	      }
	    else // slim().degree() == LINEAR
	      {
		Linear<T>& bf = (Linear<T>&) sballs[j].bf();
		std::vector<T> c( LINEAR_COEFF_ALL );
		bf.getCoeffs( c );
		sballs[j].getCenter( &c[0] );
		sballs[j].getSupport( &(c[7]) );
#if 0
		out.write( (char*) &c[0], sizeof( c ) );
#endif
#if 1
		fwrite( &c[0], sizeof(float), LINEAR_COEFF_ALL, out );
#endif
	      
	      }
	  }
      }

#if 0
    out.close();
#endif
#if 1
    fclose( out );
#endif

    return true;
  };


private:

  Slim<T>* slim_;
  
};

typedef SlimIO<double> SlimIOd;
typedef SlimIO<float>  SlimIOf;

#endif // _BALLIO_HXX

