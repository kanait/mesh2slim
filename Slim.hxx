////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIM_HXX
#define _SLIM_HXX 1

#include <vector>
using namespace std;

#include "SlimBall.hxx"

template <typename T>
class Slim {

public:

  Slim() {
    init();
  };
  ~Slim() {
    clear();
  };

  void init() {
    oriented_ = false;
  };

  void clear() {
    for ( int i = 0; i < sbarrays_.size(); ++i ) sbarrays_[i].clear(); 
    sbarrays_.clear();
  };
  
  bool oriented() const { return oriented_; };
  void setOriented( bool f ) { oriented_ = f; };

  int degree() const { return degree_; };
  void setDegree( int d ) { degree_ = d; };

  void setBallArraySize( int n ) { sbarrays_.resize( n ); };

  int level() const { return level_; };
  void setLevel( int l ) { level_ = l; };

  std::vector< SlimBall<T> >& sbarray( int i ) { return sbarrays_[i]; };

  void print() {
    
    int bcount = 0;
    int blcount = 0;
    for ( int i = 0; i < sbarrays_.size(); ++i )
      {
	for ( int j = 0; j < sbarrays_[i].size(); ++j )
	  {
	    bcount++;
	    if ( sbarrays_[i][j].isLeaf() ) blcount++;
	  }
      }
    
    std::cout << "slim: "
	      << "degree " << degree_ 
	      << " level " << level_ 
	      << " oriented " << oriented_ 
	      << " balls " << bcount 
	      << " leaf " << blcount 
	      << std::endl;
  };

private:

  int degree_;
  int level_;
  bool oriented_;

  std::vector< std::vector< SlimBall<T> > > sbarrays_;

};

typedef Slim<double> Slimd;
typedef Slim<float>  Slimf;

#endif // _SLIM_HXX
