////////////////////////////////////////////////////////////////////
//
// $Id: MeshR.hxx 2021/06/05 12:38:07 kanai Exp $
//
// Mesh class for Rendering (double)
//
// Copyright (c) 2021 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _MESHR_HXX
#define _MESHR_HXX 1

#include "envDep.h"
#include "mydef.h"

#include <vector>
#include <map>
using namespace std;

#include "Point3.h"
#include "Vector3.h"
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#define ASSIGN_VERTEX 0
#define ASSIGN_FACE   1

// #define CREASE_ANGLE (M_PI / 18.0)
#define CREASE_ANGLE (M_PI / 2.0)

template <typename T>
class MeshR {

public:
  
  MeshR() { clear(); };
  ~MeshR() {};

  void clear() {

    points_.clear();
    normals_.clear();
    texcoords_.clear();
    colors_.clear();
    indices_.clear();
    fnormals_.clear();
    fareas_.clear();
    face_mates_.clear();

    n_points_ = 0;
    n_normals_ = 0;
    n_texcoords_ = 0;
    n_colors_ = 0;
    n_indices_ = 0;
    n_tex_ = 2; // x, y coordinates
  };

  // get real number of elements
  int numPoints() const { return (int) ((float) points_size() / (float) nXYZ); };
  int numNormals() const { return (int) ((float) normals_size() / (float) nXYZ); };
  int numTexcoords() const { return (int) ((float) texcoords_size() / (float) n_tex_); };
  int numColors() const { return (int) ((float) colors_size() / (float) nXYZ); };
  int numFaces() const { return (int) ((float) indices_size() / (float) TRIANGLE); };

  // get the number of internal elements
  int points_size() const { return n_points_; };
  int normals_size() const { return n_normals_; };
  int texcoords_size() const { return n_texcoords_; };
  int colors_size() const { return n_colors_; };
  int indices_size() const { return n_indices_; };

  std::vector<T>& points() { return points_; };
  std::vector<T>& normals() { return normals_; };
  std::vector<T>& texcoords() { return texcoords_; };
  std::vector<T>& colors() { return colors_; };
  std::vector<unsigned int>& color_ids() { return color_ids_; };
  std::vector<unsigned int>& indices() { return indices_; };
  std::vector<T>& fnormals() { return fnormals_; };
  std::vector<T>& fareas() { return fareas_; };
  std::vector<int>& face_mates() { return face_mates_; };

  Point3<T>& point( int i ) { return Point3<T>( points_[3*i], points_[3*i+1], points_[3*i+2] ); };
  Vector3<T>& normal( int i ) { return Vector3<T>( normals_[3*i], normals_[3*i+1], normals_[3*i+2] ); };
  Point3<T>& texcoord( int i ) { return Point3<T>( texcoords_[3*i], texcoords_[3*i+1], texcoords_[3*i+2] ); };
  Point3<T>& color( int i ) { return Point3<T>( colors_[3*i], colors_[3*i+1], colors_[3*i+2] ); };
  int index( int i ) const { return indices_[i]; };
  T fnormal( int i ) const { return fnormals_[i]; };
  void fnormal( int i, Vector3<T>& nrm ) const { nrm.set( fnormals_[nXYZ*i],
							  fnormals_[nXYZ*i+1],
							  fnormals_[nXYZ*i+2] ); };
  T farea( int i ) const { return fareas_[i]; };
  int face_mate( int i ) const { return face_mates_[i]; };

  void clearFaceNormals() { fnormals_.clear(); };
  void clearFaceAreas() { fareas_.clear(); };
  void clearFaceMates() { face_mates_.clear(); };

  int n_tex() const { return n_tex_; };
  void setNTex( int n ) { n_tex_ = n; };

  void setColorAssigned( unsigned short n ) { color_assigned_ = n; };
  unsigned short colorAssigned() const { return color_assigned_; };

  //
  // "add" elements
  //

  void addIndex( unsigned int i ) {
    indices_.push_back( i ); n_indices_++;
  };

  void addPoint( T x, T y, T z ) {
    points_.push_back( x ); n_points_++;
    points_.push_back( y ); n_points_++;
    points_.push_back( z ); n_points_++;
  };

  void addNormal( T x, T y, T z ) {
    normals_.push_back( x ); n_normals_++;
    normals_.push_back( y ); n_normals_++;
    normals_.push_back( z ); n_normals_++;
  };

  void addTexcoord( T x, T y ) {
    n_tex_ = 2;
    texcoords_.push_back( x ); n_texcoords_++;
    texcoords_.push_back( y ); n_texcoords_++;
  };

  void addTexcoord( T x, T y, T z ) {
    n_tex_ = 3;
    texcoords_.push_back( x ); n_texcoords_++;
    texcoords_.push_back( y ); n_texcoords_++;
    texcoords_.push_back( z ); n_texcoords_++;
  };

  void addColor( T x, T y, T z ) {
    colors_.push_back( x ); n_colors_++;
    colors_.push_back( y ); n_colors_++;
    colors_.push_back( z ); n_colors_++;
  };
  
  //
  // "set" elements
  //

  void setPoint( int i, Point3<T>& p ) {
    points_[i]   = p.x;
    points_[i+1] = p.y;
    points_[i+2] = p.z;
  };

  void setPoint( int i, T x, T y, T z ) {
    points_[i]   = x;
    points_[i+1] = y;
    points_[i+2] = z;
  };

  void setNormal( int i, Vector3<T>& p ) {
    normals_[i]   = p.x;
    normals_[i+1] = p.y;
    normals_[i+2] = p.z;
  };

  void setNormal( int i, T x, T y, T z ) {
    normals_[i]   = x;
    normals_[i+1] = y;
    normals_[i+2] = z;
  };
    
  void setTexcoord( int i, Point3<T>& p ) {
    texcoords_[i]   = p.x;
    texcoords_[i+1] = p.y;
    texcoords_[i+2] = p.z;
  };

  void setTexcoord( int i, T x, T y, T z ) {
    texcoords_[i]   = x;
    texcoords_[i+1] = y;
    texcoords_[i+2] = z;
  };

  void setTexcoord( int i, T x, T y ) {
    texcoords_[i]   = x;
    texcoords_[i+1] = y;
  };

  void setColor( int i, Point3<T>& p ) {
    colors_[i]   = p.x;
    colors_[i+1] = p.y;
    colors_[i+2] = p.z;
  };

  void setColor( int i, T x, T y, T z ) {
    colors_[i]   = x;
    colors_[i+1] = y;
    colors_[i+2] = z;
  };

  void setColorId( int i, unsigned int n ) {
    color_ids_[i] = n;
  };

  void setIndex( int i, unsigned int f ) {
    indices_[i] = f;
  };

  void setFaceMate( int i, unsigned int f ) {
    face_mates_[i] = f;
  };

  void deleteTexcoords() { texcoords_.clear(); n_texcoords_ = 0; };
  
  void reservePoints( int n )  { n_points_  = n * nXYZ;     points_.resize( n_points_ ); };
  void reserveNormals( int n ) { n_normals_ = n * nXYZ;     normals_.resize( n_normals_ ); };
  void reserveTexcoords( int n, int t ) { n_texcoords_ = n * t; texcoords_.resize( n_texcoords_ ); };
  void reserveColors( int n ) { n_colors_ = n * nXYZ;       colors_.resize( n_colors_ ); };
  void reserveColorIds( int n ) { n_color_ids_ = n;         color_ids_.resize( n_color_ids_ ); };
  void reserveIndices( int n ) { n_indices_ = n * TRIANGLE; indices_.resize( n_indices_ ); };
  void reserveFaceMates( int n ) { face_mates_.resize( n * TRIANGLE ); };

  void getFacePoints( unsigned int face_id, Point3<T>& p0, Point3<T>& p1, Point3<T>& p2 ) {
    unsigned int i0 = nXYZ * indices_[ TRIANGLE * face_id ];
    unsigned int i1 = nXYZ * indices_[ TRIANGLE * face_id + 1 ];
    unsigned int i2 = nXYZ * indices_[ TRIANGLE * face_id + 2 ];
    p0.set( points_[ i0 ], points_[ i0 + 1 ], points_[ i0 + 2 ] );
    p1.set( points_[ i1 ], points_[ i1 + 1 ], points_[ i1 + 2 ] );
    p2.set( points_[ i2 ], points_[ i2 + 1 ], points_[ i2 + 2 ] );
  };

  void getFaceMates( unsigned int face_id, unsigned int m[TRIANGLE] ) {
    m[0] = (int) ( (float) face_mates_[ TRIANGLE * face_id ] / 3.0f );
    m[1] = (int) ( (float) face_mates_[ TRIANGLE * face_id + 1 ] / 3.0f );
    m[2] = (int) ( (float) face_mates_[ TRIANGLE * face_id + 2 ] / 3.0f );
  };

  // 縮退チェック機能付
  bool checkFaceNormal( const unsigned int face_id ) {

    Point3<T> p0;
    Point3<T> p1;
    Point3<T> p2;
    getFacePoints( face_id, p0, p1, p2 );

    Vector3<T> v1( p1 - p0 );
    Vector3<T> v2( p2 - p0 );
    Vector3<T> nrm;
    nrm.cross( v1, v2 );
    if ( nrm.length() < 1.0e-07 )
      {
	return false;
      }

    return true;
  };

  void calcNormal( Point3<T>& p0, Point3<T>& p1, Point3<T>& p2, Vector3<T>& nrm,
		   const bool isNormalized = true ) {
    Vector3<T> v1( p1 - p0 );
    Vector3<T> v2( p2 - p0 );
    nrm.cross( v1, v2 );
    if ( isNormalized == true )
      nrm.normalize();
  };

  T calcFaceArea( Point3<T>& p0, Point3<T>& p1, Point3<T>& p2 ) {
    Vector3<T> v1( p1 - p0 );
    Vector3<T> v2( p2 - p0 );
    Vector3<T> a;
    a.cross(v1, v2);
    return a.length() * .5f;
  };

  T dihedralAngle( int id0, int id1 ) {
    Vector3<T> n0, n1;
    fnormal( id0, n0 );
    fnormal( id1, n1 );
    return n0.angle( n1 );
  };


  void scale( T f ) {
    for ( int i = 0; i < points_.size(); i += nXYZ )
      {
	points_[i] *= f;
	points_[i+1] *= f;
	points_[i+2] *= f;
      }
  };

  void normalize() {
    std::cout << "normalize ... " << endl;
    std::cout << "points size = " << points_.size() << endl;

    Point3<T> vmax, vmin;
//     int i = 0;
    for ( int i = 0; i < points_.size(); i += nXYZ )
      {
	T x = points_[i]; 
	T y = points_[i+1]; 
	T z = points_[i+2];
	if ( i )
	  {
	    if (x > vmax.x) vmax.x = x;
	    if (x < vmin.x) vmin.x = x;
	    if (y > vmax.y) vmax.y = y;
	    if (y < vmin.y) vmin.y = y;
	    if (z > vmax.z) vmax.z = z;
	    if (z < vmin.z) vmin.z = z;
	  }
	else
	  {
	    vmax.set( x, y, z ); vmin.set( x, y, z );
	  }
      }

    Point3<T> cen = vmax + vmin; cen.scale(.5);
    Point3<T> len = vmax - vmin;
    T maxlen = (std::fabs(len.x) > std::fabs(len.y) )
      ? std::fabs(len.x) : std::fabs(len.y); 
    maxlen = ( maxlen > std::fabs(len.z) ) ? maxlen : std::fabs(len.z);
    
    for ( int i = 0; i < points_.size(); i += nXYZ )
      {
	T x = (points_[i]   - cen.x) / maxlen;
	T y = (points_[i+1] - cen.y) / maxlen;
	T z = (points_[i+2] - cen.z) / maxlen;
	points_[i]   = x;
	points_[i+1] = y;
	points_[i+2] = z;
      }

    std::cout << "done." << endl;
    
  };
  
  void createFaceNormals( const bool isNormalized = true ) {
    std::cout << "calc face normals ... " << endl;
    fnormals_.resize( numFaces() * nXYZ );
    for ( int i = 0; i < numFaces(); ++i )
      {
	Vector3<T> nrm;
	Point3<T> p0( points_[ nXYZ * indices_[ TRIANGLE * i ] ],
		      points_[ nXYZ * indices_[ TRIANGLE * i ] + 1 ],
		      points_[ nXYZ * indices_[ TRIANGLE * i ] + 2 ] );
	Point3<T> p1( points_[ nXYZ * indices_[ TRIANGLE * i + 1 ] ],
		      points_[ nXYZ * indices_[ TRIANGLE * i + 1 ] + 1 ],
		      points_[ nXYZ * indices_[ TRIANGLE * i + 1 ] + 2 ] );
	Point3<T> p2( points_[ nXYZ * indices_[ TRIANGLE * i + 2 ] ],
		      points_[ nXYZ * indices_[ TRIANGLE * i + 2 ] + 1 ],
		      points_[ nXYZ * indices_[ TRIANGLE * i + 2 ] + 2 ] );
	calcNormal( p0, p1, p2, nrm, isNormalized );
	fnormals_[ nXYZ * i ] = nrm.x;
	fnormals_[ nXYZ * i + 1 ] = nrm.y;
	fnormals_[ nXYZ * i + 2 ] = nrm.z;
      }
    std::cout << "done."<< std::endl;
  };

  void createFaceAreas() {
    std::cout << "calc face areas ... " << endl;
    fareas_.resize( numFaces() );
    for ( int i = 0; i < numFaces(); ++i )
      {
	unsigned int i0 = indices_[ TRIANGLE * i ];
	unsigned int i1 = indices_[ TRIANGLE * i + 1 ];
	unsigned int i2 = indices_[ TRIANGLE * i + 2 ];
	Point3<T> p0( points_[ nXYZ * i0 ],
		      points_[ nXYZ * i0 + 1 ],
		      points_[ nXYZ * i0 + 2 ] );
	Point3<T> p1( points_[ nXYZ * i1 ],
		      points_[ nXYZ * i1 + 1 ],
		      points_[ nXYZ * i1 + 2 ] );
	Point3<T> p2( points_[ nXYZ * i2 ],
		      points_[ nXYZ * i2 + 1 ],
		      points_[ nXYZ * i2 + 2 ] );
	fareas_[ i ] = calcFaceArea( p0, p1, p2 );
      }
    std::cout << "done."<< std::endl;
  };

  void createVertexNormals() {
    std::cout << "calculate normals ... " << std::endl;

    //   std::cout << "points " << points_.size() << std::endl;

    normals_.resize( points_.size() );
  
    for ( int i = 0; i < normals_.size(); ++i ) normals_[i] = .0f;
    std::vector<T> weights; weights.resize(  numPoints() );
    for ( int i = 0; i < weights.size(); ++i ) weights[i] = .0f;

    // store face normals multiplied by area weights
    for ( int i = 0; i < indices_.size(); i += TRIANGLE )
      {
	unsigned int id0 = indices_[i];
	unsigned int id1 = indices_[i+1];
	unsigned int id2 = indices_[i+2];
	Point3<T> p0 ( points_[ nXYZ * id0 ], 
		       points_[ nXYZ * id0 + 1 ],
		       points_[ nXYZ * id0 + 2 ] );
	Point3<T> p1 ( points_[ nXYZ * id1 ], 
		       points_[ nXYZ * id1 + 1 ],
		       points_[ nXYZ * id1 + 2 ] );
	Point3<T> p2 ( points_[ nXYZ * id2 ], 
		       points_[ nXYZ * id2 + 1 ],
		       points_[ nXYZ * id2 + 2 ] );
      
	Vector3<T> nrm;
	calcNormal( p0, p1, p2, nrm );
	T w = calcFaceArea( p0, p1, p2 );
	nrm.scale( w );

	normals_[ nXYZ * (id0) ]     += nrm.x;
	normals_[ nXYZ * (id0) + 1 ] += nrm.y;
	normals_[ nXYZ * (id0) + 2 ] += nrm.z;
	weights[ id0 ] += w;
      
	normals_[ nXYZ * (id1) ]     += nrm.x;
	normals_[ nXYZ * (id1) + 1 ] += nrm.y;
	normals_[ nXYZ * (id1) + 2 ] += nrm.z;
	weights[ id1 ] += w;

	normals_[ nXYZ * (id2) ]     += nrm.x;
	normals_[ nXYZ * (id2) + 1 ] += nrm.y;
	normals_[ nXYZ * (id2) + 2 ] += nrm.z;
	weights[ id2 ] += w;
      }

    // divided by weight and normalize
    for ( int i = 0; i < numPoints(); ++i )
      {
	normals_[ nXYZ * i ]     /= weights[i];
	normals_[ nXYZ * i + 1 ] /= weights[i];
	normals_[ nXYZ * i + 2 ] /= weights[i];

	Vector3<T> nrm( normals_[ nXYZ * i ], normals_[ nXYZ * i + 1 ], normals_[ nXYZ * i + 2 ] );
	nrm.normalize();

	normals_[ nXYZ * i ]     = nrm.x;
	normals_[ nXYZ * i + 1 ] = nrm.y;
	normals_[ nXYZ * i + 2 ] = nrm.z;
      }

    n_normals_ = normals_.size();

    std::cout << "done. " << std::endl;
    
  };

  void createFaceMates() {
    std::cout << "create index pair map ... " << std::endl;

    typedef std::pair< unsigned int, unsigned int > type_uu;
    typedef std::pair< unsigned int, std::pair< unsigned int, unsigned int> > type_up;
    typedef std::multimap< unsigned int, std::pair< unsigned int, unsigned int> >::iterator type_upIter;
    // multimap< id0, pair< id1, fid > >
    std::multimap< unsigned int, std::pair< unsigned int, unsigned int> > indices_pair;

    // 領域の確保
    face_mates_.resize( indices_.size() );

    for( int i = 0; i < indices_.size(); i += TRIANGLE )
      {
#if 1
	unsigned int id[3];
	id[0] = indices_[ i ];
	id[1] = indices_[ i + 1 ];
	id[2] = indices_[ i + 2 ];

	indices_pair.insert( type_up(id[0], type_uu(id[1], i)) );
	indices_pair.insert( type_up(id[1], type_uu(id[2], i)) );
	indices_pair.insert( type_up(id[2], type_uu(id[0], i)) );
#endif

#if 0
	unsigned int id[3];
	id[0] = indices_[ i ];
	id[1] = indices_[ i + 1 ];
	id[2] = indices_[ i + 2 ];

	type_upIter cIterI[3];
	cIterI[0] = indices_pair.insert( type_up(id[0], type_uu(id[1], i)) );
	cIterI[1] = indices_pair.insert( type_up(id[1], type_uu(id[2], i)) );
	cIterI[2] = indices_pair.insert( type_up(id[2], type_uu(id[0], i)) );

	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    unsigned int id0 = indices_[ i + j ];
	    unsigned int id1 = ( j != 2 ) ? indices_[ i + (j + 1) ] : indices_[ i ];

	    // id1 をキーとする pair を見つける
	    std::pair< type_upIter, type_upIter > cIterPair
	      = indices_pair.equal_range( id1 );
	  
	    face_mates_[ i + j ] = -1;
	    bool found = false;
	    for ( type_upIter cIter = cIterPair.first; cIter != cIterPair.second; ++cIter )
	      {
		// pair のうち 最初のキーが id0 のものを見つける
		if ( (*cIter).second.first == id0 ) // found!
		  {
		    found = true;
		    face_mates_[ i + j ] = (*cIter).second.second;
		    face_mates_[ (*cIter).second.second ] = i;

		    indices_pair.erase( cIter );
		    indices_pair.erase( cIterI[j] );
		    break;
		  }
	      }
#if 0
	    if ( found == false )
	      {
		indices_pair.insert( type_up(id0, type_uu(id1, i)) );
	      }
#endif
	  }
#endif
      }

#if 1
    // メイトの作成
    std::cout << "create face mates ... " << std::endl;
  
    for( int i = 0; i < indices_.size(); i += TRIANGLE )
      {
	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    unsigned int id0 = indices_[ i + j ];
	    unsigned int id1 = ( j != 2 ) ? indices_[ i + (j + 1) ] : indices_[ i ];

	    // id1 をキーとする pair を見つける
	    std::pair< type_upIter, type_upIter > cIterPair
	      = indices_pair.equal_range( id1 );
	  
	    face_mates_[ i + j ] = -1;
	    for ( type_upIter cIter = cIterPair.first; cIter != cIterPair.second; ++cIter )
	      {
		// pair のうち 最初のキーが id0 のものを見つける
		if ( (*cIter).second.first == id0 ) // found!
		  {
		    face_mates_[ i + j ] = (*cIter).second.second;
		  }
	      }
	  }
      }

#endif

    std::cout << "done. " << std::endl;
    
  };

  void deleteFaceMates() { face_mates_.clear(); };

  void printInfo() {
    std::cout << "mesh " << " ";
    if ( points_.size() ) std::cout << " v " << numPoints() << " ";
    if ( normals_.size() ) std::cout << " n " << numNormals() << " ";
    if ( texcoords_.size() ) std::cout << " t " << numTexcoords() << " ";
    if ( indices_.size() ) std::cout << " f " << numFaces() << " ";
    std::cout << std::endl;
  };

  void getBBox( Point3<T>& vmax, Point3<T>& vmin ) {
    for ( int i = 0; i < points_.size(); i += 3 )
      {
	Point3<T> p( points_[i], points_[i+1], points_[i+2] );
	if ( i )
	  {
	    if (p.x > vmax.x) vmax.x = p.x;
	    if (p.x < vmin.x) vmin.x = p.x;
	    if (p.y > vmax.y) vmax.y = p.y;
	    if (p.y < vmin.y) vmin.y = p.y;
	    if (p.z > vmax.z) vmax.z = p.z;
	    if (p.z < vmin.z) vmin.z = p.z;
	  }
	else
	  {
	    vmax.set( p ); vmin.set( p );
	  }
      }
  };


  
private:

  // vertex points
  int n_points_;
  std::vector<T> points_;

  // normals
  int n_normals_;
  std::vector<T> normals_;
 
  // texture coordinates
  int n_tex_;
  int n_texcoords_;
  std::vector<T> texcoords_;

  // colors
  unsigned short color_assigned_;
  int n_colors_;
  std::vector<T> colors_;

  // assigned color for vertices or faces
  int n_color_ids_;
  std::vector<unsigned int> color_ids_;

  // face indexes
  int n_indices_;
  std::vector<unsigned int> indices_;

//   int n_indices_;
//   std::vector<unsigned int> normal_indices_;

  // face normals
  std::vector<T> fnormals_;

  // face areas
  std::vector<T> fareas_;

  // face mates
  std::vector<int> face_mates_;

};

typedef MeshR<double> MeshRd;
typedef MeshR<float>  MeshRf;

#endif // _MESHR_HXX

