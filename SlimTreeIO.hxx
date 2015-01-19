////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIMTREEIO_HXX
#define _SLIMTREEIO_HXX 1

#include <vector>
#include <cstdio>
#include <iostream>
#include <fstream>
using namespace std;

#include "BasisFunction.hxx"
#include "SlimBall.hxx"
#include "SlimTree.hxx"
#include "SlimTreeApp.hxx"

template <typename T>
class SlimTreeIO : public SlimTreeApp<T> {

public:

  SlimTreeIO( SlimTree<T>& slim ) { setSlimTree( slim ); };
  ~SlimTreeIO(){};

  bool inputFromFile( const char* const infile ) {
    
    FILE* in = std::fopen( infile, "r+b" );

    int N[3];
    std::fread( N, sizeof(int), 3, in );

    bool oriented;
    int degree = N[0];
    int num_slimballs  = N[1];
    std::cout << "#slimballs: " << num_slimballs << std::endl;

    if( degree < 0 )
      {
	degree = -degree;
	oriented = false;
      }
    else
      oriented = true;
    int num_user_defined = N[2];
  
    slimtree().setDegree( degree );
    slimtree().setOriented( oriented );
    slimtree().setNumUserDefined( num_user_defined );

    T minDis = .0;
    T maxDis = .0;

    int id = 0;
    for ( int i = 0; i < num_slimballs; ++i )
      {
	SlimBall<T>* slimball = slimtree().addSlimBall();
	slimball->setID( id ); ++id;

	if ( degree == CUBIC )
	  {
#if 0
	    std::vector<T> c( CUBIC_COEFF_ALL );

	    std::fread( &c[0], sizeof(T), CUBIC_COEFF_ALL, in );

	    sballs[j].setCenter( c[0], c[1], c[2] );
	    sballs[j].setSupport( c[23] );
	      
	    Cubic<T>* bf = new Cubic<T>;
	    bf->setCoeffs( c );
	    sballs[j].addBasisFunction( *bf );
#endif
	  }
	else if ( degree == QUADRATIC )
	  {
	    int num_param = 3 + 1 + QUADRATIC_COEFF + 1 + num_user_defined;
	    
	    // 一旦 float で読み込んで，あとで T に変換する
	    std::vector<float> cc( num_param );

	    std::fread( &cc[0], sizeof(float), num_param, in );

	    // float から T に変換
	    std::vector<T> c( num_param );
	    for ( int j = 0; j < num_param; ++j ) c[j] = (T) cc[j];

	    // center
	    slimball->setCenter( c[0], c[1], c[2] );
	    
	    // 多項式の係数
	    Quadratic<T>* bf = new Quadratic<T>;
	    bf->setCoeffs( c );
	    slimball->addBasisFunction( *bf );

	    // radius
	    slimball->setSupport( c[13] );

	    // parent
	    int parent_id = (int) c[14];
	    if ( parent_id != -1 )
	      {
		SlimBall<T>* parent = slimtree().slimball( parent_id - 1 );
		slimball->setParent( parent );
		parent->addChild( slimball );
// 		slimball->setLevel( parent->level() + 1 );
	      }
	    else
	      {
// 		slimball->setLevel( 0 );
	      }
	    
	    // user defined values
	    slimball->setNumUserDefined( num_user_defined );
	    for ( int j = 0; j < num_user_defined; ++j )
	      {
		slimball->setUserDefined( j, c[15+j] );
	      }
	    
	    if ( i ) 
	      {
		if ( slimball->userDefined() < minDis ) minDis = slimball->userDefined();
		if ( slimball->userDefined() > maxDis ) maxDis = slimball->userDefined();
	      }
	    else
	      {
		minDis = slimball->userDefined();
		maxDis = slimball->userDefined();
	      }

	  }
	else
	  { 
#if 0
	    std::vector<T> c( LINEAR_COEFF_ALL );

	    std::fread( &c[0], sizeof(T), LINEAR_COEFF_ALL, in );

	    sballs[j].setCenter( c[0], c[1], c[2] );
	    sballs[j].setSupport( c[7] );

	    Linear<T>* bf = new Linear<T>;
	    bf->setCoeffs( c );
	    sballs[j].addBasisFunction( *bf );
#endif
	  }
      }
    
    std::fclose( in );

    cout << "distance min: " << minDis << " max: " << maxDis << endl;

    return true;
  };

  bool outputToFile( const char* const outfile ) {

    cout << "save slimtree file: " << outfile << " ... " << endl;

    // set IDs for save
    std::vector< SlimBall<T>* > slimballs = slimtree().slimballs();
    for( int i = 0; i < slimballs.size(); ++i )
      {
	slimballs[i]->setID( i + 1 );
      }

    FILE* out = fopen( outfile, "w+b" );

    int L[3];
    L[0] = ( slimtree().oriented() ) ? slimtree().degree() : -slimtree().degree();
    L[1] = slimtree().slimballs().size();
    L[2] = slimtree().numUserDefined();

    fwrite( L, sizeof(int), 3, out );

    for( int i = 0; i < slimballs.size(); ++i )
      {
	if ( slimtree().degree() == CUBIC )
	  {
#if 0
	    Cubic<T>& bf = (Cubic<T>&) sballs[j].bf();
	    std::vector<T> c( CUBIC_COEFF_ALL );
	    bf.getCoeffs( c );
	    sballs[j].getCenter( &c[0] );
	    sballs[j].getSupport( &(c[23]) );

	    fwrite( &c[0], sizeof(float), CUBIC_COEFF_ALL, out );
#endif
	  }
	else if ( slimtree().degree() == QUADRATIC )
	  {
	    int num_param = 3 + 1 + QUADRATIC_COEFF + 1 + slimtree().numUserDefined();

	    std::vector<T> c( num_param );
	    
	    // center (0-2)
	    slimballs[i]->getCenter( &c[0] );

	    // quadratic coefficients (3-12)
	    Quadratic<T>& bf = (Quadratic<T>&) slimballs[i]->bf();
	    bf.getCoeffs( c );
	    
	    // support radius (13)
	    slimballs[i]->getSupport( &(c[13]) );

	    // parents (14)
	    c[14] = ( slimballs[i]->parent() ) ? slimballs[i]->parent()->id() : -1;

	    // ユーザ定義変数 (15-)
	    std::vector<T>& ud = slimballs[i]->userDefineds();
	    for ( int j = 0; j < ud.size(); ++j ) c[15+j] = ud[j];

	    // float へ変換
	    std::vector<float> cc( num_param );
	    for ( int k = 0; k < num_param; ++k )
	      {
		cc[k] = (float) c[k];
	      }

	    fwrite( &cc[0], sizeof(float), num_param, out );

	  }
	else // slimtree().degree() == LINEAR
	  {
#if 0
	    Linear<T>& bf = (Linear<T>&) sballs[j].bf();
	    std::vector<T> c( LINEAR_COEFF_ALL );
	    bf.getCoeffs( c );
	    sballs[j].getCenter( &c[0] );
	    sballs[j].getSupport( &(c[7]) );

	    fwrite( &c[0], sizeof(float), LINEAR_COEFF_ALL, out );
#endif	      
	  }
      }

    fclose( out );

    cout << "save slimtree file: done." << endl;
    return true;
  };

};

typedef SlimTreeIO<double> SlimTreedIO;
typedef SlimTreeIO<float>  SlimTreefIO;

#endif // _SLIMTREEIO_HXX

