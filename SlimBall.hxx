////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIMBALL_HXX
#define _SLIMBALL_HXX 1

#include <vector>
using namespace std;

#include <Point3.h>
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

//for spline patched Gaussian
#define EE16 2.1653645317858030703

#include "BasisFunction.hxx"

#define DEFAULT_NUM_USER_DEFINED 1

// for distance
#define DIST_AVR 0
#define DIST_MIN 1

template <typename T>
class SlimBall {

public:

  SlimBall() : parent_(NULL), id_(-1) { init(); };
  ~SlimBall() {
    deleteBases();
    userDefined_.clear();
  };

  void init() {
    setNumUserDefined( DEFAULT_NUM_USER_DEFINED );
  };

  void setID( int i ) { id_ = i; };
  int id() const { return id_; };

  T support() const { return support_; };
  T supportSquared() const { return support_*support_; };
  void setSupport( T f ) { support_ = f; };
  void getSupport( T* x ) { *x = support_; };

  void setCenter( T x, T y, T z ) { center_.set( x, y, z ); };
  void setCenter( Point3<T>& p ) { center_.set( p ); };
  Point3<T>& center() { return center_; };
  void getCenter( T* p ) {
    p[0] = center_.x;
    p[1] = center_.y;
    p[2] = center_.z;
  };
  void getCenter( Point3<T>& p ) {
    p.x = center_.x;
    p.y = center_.y;
    p.z = center_.z;
  };

//   void setIsLeaf( bool f ) { isLeaf_ = f; };
  bool isLeaf() const { return ( childs_.size() ) ? false : true; };

  void deleteBases() {
    for ( int i = 0; i < bases_.size(); ++i ) delete bases_[i];
    bases_.clear();
  };
  //BasisFunction<T>& bf() { return *(bases_[0]); };
  BasisFunction<T>& bf( int i = 0 ) const { return *(bases_[i]); };
  void addBasisFunction( BasisFunction<T>& b ) { bases_.push_back( &b ); };

  void setParent( SlimBall<T>* sb ) { parent_ = sb; };
  SlimBall<T>* parent() { return parent_; };

  void addChild( SlimBall<T>* sb ) { childs_.push_back( sb ); };
  std::vector< SlimBall<T>* >& childs() { return childs_; };
  SlimBall<T>* child( int i ) { return childs_[i]; };

  // user defined value functions
  void setNumUserDefined( int n ) { userDefined_.resize( n ); };
  std::vector<T>& userDefineds() { return userDefined_; };
  void setUserDefined( T f ) { userDefined_[0] = f; };
  void setUserDefined( int i, T f ) { userDefined_[i] = f; };
  T userDefined() const { return userDefined_[0]; };
  T userDefined( int i ) const { return userDefined_[i]; };
  

  T weight( T x, T y, T z ) const { 

    weight( Point3<T>( x, y, z ) );

  }

  // p は絶対座標
  T weight( Point3<T>& p ) const {

    T d = p.distanceSquared( center_ );
    if ( d > supportSquared() )
      {
	return .0f;
      }
    else
      {
	return weight_( std::sqrt(d), support_ );
      }

  };
  
  int countBalls( T error ) const {

    int count = 0;
    if ( userDefined() < error ) 
      {
	++count;
	return count;
      }

    for ( int i = 0; i < childs_.size(); ++i )
      {
	count += childs_[i]->countBalls( error );
      }

    return count;
  };

  double distance( Point3<T>& p, T error, int* n, int flag ) const {

    *n = 0;
    
    if ( center_.distance( p ) > support_ ) return -1;

    BasisFunction<T>& b = (BasisFunction<T>&) bf();

    T dis;
    if ( flag == DIST_AVR ) dis = 0;
    else dis = -1;

    if ( (childs_.empty()) || (userDefined() < error) )
      {
	Point3<T> psub( p - center_ );
	T a = std::fabs( b.poly( psub ) );
	Vector3<T> g; b.polyG( psub, g );
	T c = std::sqrt( g.x*g.x+g.y*g.y+g.z*g.z );
	dis = a / c;
	*n = 1;
// 	if (childs_.empty())
// 	cout << "leaf node c0 " << b.coeff(0) << " c1 " << b.coeff(1) << " c2 " << b.coeff(2) 
// 	     << " dis " << a << endl;
	return dis;
      }

    bool fflag = false;
    for ( int i = 0; i < childs_.size(); ++i )
      {
	int n1;
	T dis0 = childs_[i]->distance( p, error, &n1, flag );
	if ( dis0 < -0.1 ) continue;

	if ( flag == DIST_AVR )
	  {
	    dis += dis0;
	  }
	else // DIST_MIN
	  {
	    if ( !fflag ) 
	      {
		dis = dis0;
		fflag = true;
	      }
	    else
	      {
		if ( dis0 < dis ) dis = dis0;
	      }
	  }
	
	*n += n1;
      }

    return dis;
  };

  double grad_distance( Point3<T>& p, Vector3<T>& nrm, T error, int* n, int flag ) const {

    *n = 0;
    
    if ( center_.distance( p ) > support_ ) return -1;

    BasisFunction<T>& b = (BasisFunction<T>&) bf();

    T dis;
    if ( flag == DIST_AVR ) dis = 0;
    else dis = -1;

    if ( (childs_.empty()) || (userDefined() < error) )
      {

	Point3<T> psub( p - center_ );
	Vector3<T> g; b.polyG( psub, g );
	Vector3<T> sub( g - nrm );
	dis = sub.length() / g.length();

	*n = 1;
// 	if (childs_.empty())
// 	cout << "leaf node c0 " << b.coeff(0) << " c1 " << b.coeff(1) << " c2 " << b.coeff(2) 
// 	     << " dis " << a << endl;
	return dis;
      }

    bool fflag = false;
    for ( int i = 0; i < childs_.size(); ++i )
      {
	int n1;
	T dis0 = childs_[i]->grad_distance( p, nrm, error, &n1, flag );
	if ( dis0 < -0.1 ) continue;

	if ( flag == DIST_AVR )
	  {
	    dis += dis0;
	  }
	else // DIST_MIN
	  {
	    if ( !fflag ) 
	      {
		dis = dis0;
		fflag = true;
	      }
	    else
	      {
		if ( dis0 < dis ) dis = dis0;
	      }
	  }
	
	*n += n1;
      }

    return dis;
  };

  void print() {
    cout << "center " << center_ << " radius " << support_ << " error " << userDefined_[0] << endl;
    if ( !(bases_.empty()) )
      {
	cout << "coeff: ";
	BasisFunction<T>& b = (BasisFunction<T>&) bf();
	std::vector<T>& c = b.coeffs();
	for ( int i = 0; i < c.size(); ++i )
	  cout << c[i] << " ";
	cout << endl;
      }
  };
  
private:

  // id (for save)
  int id_;

  // user-defined value
  std::vector<T> userDefined_;

  std::vector< BasisFunction<T>* > bases_;

  Point3<T> center_;
  T   support_;

//   bool isLeaf_;

  SlimBall<T>* parent_;

  std::vector< SlimBall<T>* > childs_;

  T weight_( T d, T R ) const {
    
    if( R < d )
      return .0f;
    else
      {
	d /= R;
	if( d < .5f )
	  return std::exp( -8 * d * d );
	else
	  {
	    d = 1.0f - d;
	    d = d * d;
	    return EE16 * d * d;
	  }
      }
    
  };

};

typedef SlimBall<double> SlimBalld;
typedef SlimBall<float>  SlimBallf;

#endif // _SLIMBALL_HXX
