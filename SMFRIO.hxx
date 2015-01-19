////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2005 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SMFRIO_HXX
#define _SMFRIO_HXX 1

#include "envDep.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
using namespace std;

#include "timer.hxx"
#include "strutil.h"
#include "tokenizer.h"

#include "MeshR.hxx"
#include "RIO.hxx"

template <typename T>
class SMFRIO : public RIO<T> {

public:

  SMFRIO() : RIO<T>() {};
  SMFRIO( MeshR<T>& mesh ) : RIO<T>( mesh ) {};
  ~SMFRIO() {};
  
  bool inputFromFile( const char* const filename ) {
    //
    // file open
    //

    // 1st try
    std::ifstream ifs0; ifs0.open( filename );
    if ( !ifs0.is_open() )
      {
	std::cerr << "Cannot open " << filename << std::endl;
	return false;
      }

    std::cout << "Input from file: " << filename << endl;


    // count elements
    int v_count = 0;  // vertex
    int n_count = 0;  // normal
    int r_count = 0;  // texcoord
    int c_count = 0;  // color
    int cf_count = 0;  // color id
    int f_count = 0;  // face indices
    int fm_count = 0;  // face mates

    std::cout << "count elements ..." << std::endl;

    std::string cline;
    StrUtil strutil;
    while ( getline(ifs0, cline, '\n') )
      {
	// comment
	if ( mycomment(cline[0]) ) continue;
	else if ( cline[0] == 'v' ) ++v_count;
	else if ( cline[0] == 'n' ) ++n_count;
	else if ( cline[0] == 'r' ) 
	  {
	    if ( !r_count )
	      {
		int wc = strutil.word_count( cline );
		if ( wc == 3 ) mesh().setNTex( 2 );
		else mesh().setNTex( 3 );
	      }
	    ++r_count;
	  }
	else if ( cline[0] == 'c' ) 
	  {
	    if ( cline[1] == 'f' ) ++cf_count;
	    else ++c_count;
	  }
	else if ( cline[0] == 'f' )
	  {
	    if ( cline[1] == 'm' ) 
	      ++fm_count;
	    else 
	      ++f_count;
	  }
      }

    ifs0.close();

    std::cout << "done." << endl;

    std::cout << "meshR:";
    if ( v_count )  std::cout << " v "  << v_count;
    if ( n_count )  std::cout << " n "  << n_count;
    if ( r_count )  std::cout << " r "  << r_count;
    if ( c_count )  std::cout << " c "  << c_count;
    if ( cf_count ) std::cout << " cf " << cf_count;
    if ( f_count )  std::cout << " f "  << f_count;
    if ( fm_count )  std::cout << " fm "  << fm_count;
    std::cout << std::endl;

    std::cout << "read data .. " << std::endl;

    if ( v_count )  mesh().reservePoints(    v_count );
    if ( n_count )  mesh().reserveNormals(   n_count );
    if ( r_count )  mesh().reserveTexcoords( r_count, mesh().n_tex() );
    if ( c_count )  mesh().reserveColors(    c_count );
    if ( cf_count ) mesh().reserveColorIds(  cf_count );
    if ( f_count )  mesh().reserveIndices(   f_count );
    if ( fm_count )  mesh().reserveFaceMates(   fm_count );

  
    // 2nd try
    std::ifstream ifs; ifs.open( filename );
    if ( !ifs.is_open() )
      {
	std::cerr << "Cannot open " << filename << std::endl;
	return false;
      }

    //
    // parse smf
    //
    int v_id = 0;
    int n_id = 0;
    int r_id = 0;
    int c_id = 0;
    int cf_id = 0;
    int f_id = 0;
    int fm_id = 0;
    Timer t;
    double time0 = t.get_seconds();
    while ( getline(ifs, cline, '\n') )
      {
	std::string fw;
	strutil.first_word( cline, fw );
      
	// comment
	if ( mycomment(cline[0]) ) continue;

	// read vertices
	else if ( !fw.compare("v") )
	  {
#if 0
	    std::istringstream isstr( cline );

	    // "v"
	    std::string str; 
	    isstr >> str;
	    T x, y, z;
	    isstr >> x; isstr >> y; isstr >> z;
#endif

#if 1
	    char dummy[BUFSIZ];
	    char xc[BUFSIZ], yc[BUFSIZ], zc[BUFSIZ];
	    sscanf( cline.c_str(), "%s%s%s%s", &dummy, &xc, &yc, &zc );
	    T x = atof( xc );
	    T y = atof( yc );
	    T z = atof( zc );
#endif

	    mesh().setPoint( v_id, x, y, z );
	    v_id += nXYZ;
	  }

	// read normals
	else if ( !fw.compare("n") )
	  {
#if 0
	    std::istringstream isstr( cline );
	  
	    // "n"
	    std::string str; isstr >> str;

	    T x, y, z;
	    isstr >> x; isstr >> y; isstr >> z;
#endif

#if 1
	    char dummy[BUFSIZ];
	    char xc[BUFSIZ], yc[BUFSIZ], zc[BUFSIZ];
	    sscanf( cline.c_str(), "%s%s%s%s", &dummy, &xc, &yc, &zc );
	    T x = atof( xc );
	    T y = atof( yc );
	    T z = atof( zc );
#endif

	    mesh().setNormal( n_id, x, y, z );
	    n_id += nXYZ;
	  }

	// read texcoords
	else if ( !fw.compare("r") )
	  {
	    char dummy[BUFSIZ];
	    char xc[BUFSIZ], yc[BUFSIZ], zc[BUFSIZ];
	    T x, y, z;
	    if ( mesh().n_tex() == 2 )
	      {
		sscanf( cline.c_str(), "%s%s%s", &dummy, &xc, &yc );
		x = atof( xc );
		y = atof( yc );
		mesh().setTexcoord( r_id, x, y );
	      }
	    else if ( mesh().n_tex() == 3 )
	      {
		sscanf( cline.c_str(), "%s%s%s%s", &dummy, &xc, &yc, &zc );
		x = atof( xc );
		y = atof( yc );
		z = atof( zc );
		mesh().setTexcoord( r_id, x, y, z );
	      }
	    r_id += mesh().n_tex();
	  }

	// binding types
	else if ( !fw.compare("bind") )
	  {
	    std::string s;
	    strutil.nth_word( cline, 3, s );
	    if ( !s.compare("vertex") )
	      mesh().setColorAssigned( ASSIGN_VERTEX );
	    else if ( !s.compare("face") )
	      mesh().setColorAssigned( ASSIGN_FACE );
	  }

	// colors
	else if ( !fw.compare("c") )
	  {
	    char dummy[BUFSIZ];
	    char xc[BUFSIZ], yc[BUFSIZ], zc[BUFSIZ];
	    sscanf( cline.c_str(), "%s%s%s%s", &dummy, &xc, &yc, &zc );
	    T x = atof( xc );
	    T y = atof( yc );
	    T z = atof( zc );

	    mesh().setColor( c_id, x, y, z );
	    c_id += nXYZ;
	  }

	// face color ids
	else if ( !fw.compare("cf") )
	  {
	    char dummy[BUFSIZ];
	    char ids[BUFSIZ];
	    sscanf( cline.c_str(), "%s%s", &dummy, &ids );
	    unsigned int id = atoi( ids );

	    mesh().setColorId( cf_id, id );
	    ++ cf_id;
	  }
      
	// read faces
	else if ( !fw.compare("f") )
	  {
	    if ( cline.find("/") != std::string::npos )
	      {
		std::istringstream isstr( cline );

		// "f"
		std::string str; isstr >> str;

		int j = 0;
		while ( isstr >> str )
		  {
		    if ( j >= 3 ) break;
		    tokenizer tok( str, "/" );

		    // vertex
		    std::string val = tok.next();

		    int id0;
		    unsigned int id;
		    std::istringstream ai(val); ai >> id0;
		    id = id0 - 1;
		    mesh().setIndex( f_id, id );
	      
		    ++f_id;
	      
		    // normal
		    if ( n_count ) val = tok.next();
	      
		    ++j;
		  }
	      }
	    else
	      {
		char dummy[BUFSIZ];
		char xc[BUFSIZ], yc[BUFSIZ], zc[BUFSIZ];
		sscanf( cline.c_str(), "%s%s%s%s", &dummy, &xc, &yc, &zc );
		unsigned int id0 = atoi( xc ) - 1;
		unsigned int id1 = atoi( yc ) - 1;
		unsigned int id2 = atoi( zc ) - 1;
		mesh().setIndex( f_id, id0 ); ++f_id;
		mesh().setIndex( f_id, id1 ); ++f_id;
		mesh().setIndex( f_id, id2 ); ++f_id;
	      }
	  }
	else if ( !fw.compare("fm") )
	  {
	    char dummy[BUFSIZ];
	    char xc[BUFSIZ], yc[BUFSIZ], zc[BUFSIZ];
	    sscanf( cline.c_str(), "%s%s%s%s", &dummy, &xc, &yc, &zc );
	    unsigned int id0 = atoi( xc ) - 1;
	    unsigned int id1 = atoi( yc ) - 1;
	    unsigned int id2 = atoi( zc ) - 1;
// 	    cout << id0 << " " << id1 << " " << id2 << endl;
	    mesh().setFaceMate( fm_id, id0 * TRIANGLE ); ++fm_id;
	    mesh().setFaceMate( fm_id, id1 * TRIANGLE ); ++fm_id;
	    mesh().setFaceMate( fm_id, id2 * TRIANGLE ); ++fm_id;
	  }
      
      } // while
  
    double time1 = t.get_seconds();
    std::cout << "done. ellapsed time: " << time1 - time0 << " sec. " << std::endl;

    ifs.close();

    mesh().printInfo();
  
    return true;
  };
  
  bool outputToFile( const char* const filename,
		     bool isSaveNormal,
		     bool isSaveTexcoord, bool isSaveBLoop ) {
    setSaveNormal( isSaveNormal );
    setSaveTexcoord( isSaveTexcoord );
    setSaveBLoop( isSaveBLoop );
    return outputToFile( filename );
  };
  bool outputToFile( const char* const filename ) {
    std::ofstream ofs( filename ); if ( !ofs ) return false;

    std::cout << "Output to file: " << filename << endl;

    mesh().printInfo();

    // header
  
    ofs << "#$SMF 1.0" << std::endl;

    int vn = (int) ((float) mesh().points().size() / (float) nXYZ);
    ofs << "#$vertices " << vn << std::endl;

    int nn = (int) ((float) mesh().normals().size() / (float) nXYZ);
    if ( nn && isSaveNormal() )
      ofs << "#$normals " << nn << std::endl;

    int tn = (int) ((float) mesh().texcoords().size() / (float) mesh().n_tex());
    if ( tn && isSaveTexcoord() )
      ofs << "#$texcoords " << tn << std::endl;

    int cn  = (int) ((float) mesh().colors().size() / (float) nXYZ);
    int cin = mesh().color_ids().size();
    if ( cn && isSaveColor() )
      {
	ofs << "#$colors " << cn << std::endl;
	if ( cin ) ofs << "#$color_ids " << cin << std::endl;
      }

    int fn = (int) ((float) mesh().indices().size() / (float) TRIANGLE);
    ofs << "#$faces " << fn << std::endl;
    ofs << "#" << std::endl;

    if ( vn )
      {
	std::vector<T>& points = mesh().points();
	for ( int i = 0; i < vn; ++i )
	  {
	    ofs << "v\t"
		<< points[ nXYZ * i ] << " "
		<< points[ nXYZ * i + 1 ]  << " "
		<< points[ nXYZ * i + 2 ]  << std::endl;
	  }
      }

    if ( nn && isSaveNormal() )
      {
	std::vector<T>& normals = mesh().normals();
	for ( int i = 0; i < nn; ++i )
	  {
	    ofs << "n\t"
		<< normals[ nXYZ * i ] << " "
		<< normals[ nXYZ * i + 1 ] << " "
		<< normals[ nXYZ * i + 2 ] << std::endl;
	  }
      }

    if ( tn && isSaveTexcoord() )
      {
	int n_tex = mesh().n_tex();
	std::vector<T>& texcoords = mesh().texcoords();
	for ( int i = 0; i < tn; ++i )
	  {
	    ofs << "r\t"
		<< texcoords[ n_tex * i ] << " "
		<< texcoords[ n_tex * i + 1 ];
	    if ( n_tex == 3 )
	      {
		ofs  << " " << texcoords[ n_tex * i + 2 ];
	      }
	    ofs << std::endl;
	  }
      }

    if ( cn && isSaveColor() )
      {
	if ( mesh().colorAssigned() == ASSIGN_VERTEX )
	  {
	    ofs << "bind c vertex" << std::endl;
	  }
	else if ( mesh().colorAssigned() == ASSIGN_FACE )
	  {
	    ofs << "bind c face" << std::endl;
	  }

	// color
	std::vector<T>& colors = mesh().colors();
	for ( int i = 0; i < cn; ++i )
	  {
	    ofs << "c\t"
		<< colors[ nXYZ * i ] << " "
		<< colors[ nXYZ * i + 1 ] << " "
		<< colors[ nXYZ * i + 2 ] << std::endl;
	  }

	// color id
	std::vector<unsigned int>& color_ids = mesh().color_ids();
	for ( int i = 0; i < cin; ++i )
	  {
	    ofs << "cf\t" << color_ids[i] << endl;
	  }
      }

    if ( fn )
      {
	std::vector<unsigned int>& indices = mesh().indices();
	for ( int i = 0; i < fn; ++i )
	  {
	    ofs << "f\t" ;
	    for ( int j = 0; j < TRIANGLE; ++j )
	      {
		unsigned int id = indices[ TRIANGLE * i + j ] + 1;
		ofs << id;
		if ( nn && isSaveNormal() )   ofs << "/" << id;
		if ( tn && isSaveTexcoord() ) ofs << "/" << id;
		ofs << " " ;
	      }
	    ofs << std::endl;
	  }
      }

    ofs.close();

    return true;
    
  };

};

typedef SMFRIO<double> SMFRdIO;
typedef SMFRIO<float>  SMFRfIO;

#endif // _SMFRIO_H
