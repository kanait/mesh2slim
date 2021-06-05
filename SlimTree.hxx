////////////////////////////////////////////////////////////////////
//
// $Id: SlimTree.hxx 2021/06/05 12:42:26 kanai Exp $
//
// Copyright (c) 2021 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIMTREE_HXX
#define _SLIMTREE_HXX 1

#include "SlimBall.hxx"

#define SLIMTREE_DEGREE_DEFAULT 2
#define SLIMTREE_ORIENTED_DEFAULT false

template <typename T>
class SlimTree {

public:

  SlimTree() : degree_(SLIMTREE_DEGREE_DEFAULT), oriented_(SLIMTREE_ORIENTED_DEFAULT),
	       num_user_defined_(DEFAULT_NUM_USER_DEFINED) {};
  ~SlimTree(){
    deleteSlimBalls();
  };

  std::vector< SlimBall<T>* >& slimballs() { return slimballs_; };

  bool oriented() const { return oriented_; };
  void setOriented( bool f ) { oriented_ = f; };

  int degree() const { return degree_; };
  void setDegree( int d ) { degree_ = d; };
  
  SlimBall<T>* root() const {
    if ( !( slimballs_.empty() ) ) return slimballs_[0]; 
    else return NULL;
  };

  int numUserDefined() const { return num_user_defined_; };
  void setNumUserDefined( int i ) { num_user_defined_ = i; };

  SlimBall<T>* addSlimBall() {
    SlimBall<T>* sb = new SlimBall<T>;
    slimballs_.push_back( sb );
    return sb;
  };

  void setSlimBall( SlimBall<T>* sb ) { slimballs_.push_back( sb ); };

  SlimBall<T>* slimball( unsigned int i ) { 
    if ( i >= slimballs_.size() ) return NULL;
    return slimballs_[i]; 
  };

  void deleteSlimBalls() {
    for ( int i = 0; i < slimballs_.size(); ++i )
      {
	delete slimballs_[i];
      }
    slimballs_.clear();
  };

  double countBalls( T error ) const {
    
    if ( root() == NULL ) return -1;
    
    return root()->countBalls( error );
  };

  double distance( Point3<T>& p, T error, int flag ) const {

    if ( root() == NULL ) return -1;
    
    int n;
    double dis = root()->distance( p, error, &n, flag );

    if (!n) return -1;

    if ( flag == DIST_AVR )
      dis /= (double) n;

    return dis;
  };
  
  double grad_distance( Point3<T>& p, Vector3<T>& n, T error, int flag ) const {

    if ( root() == NULL ) return -1;
    
    int num;
    double dis = root()->grad_distance( p, n, error, &num, flag );

    if (!num) return -1;

    if ( flag == DIST_AVR )
      dis /= (double) num;

    return dis;
  };
  
private:

  int degree_;

  bool oriented_;

  int num_user_defined_;

  std::vector< SlimBall<T>* > slimballs_;

};

typedef SlimTree<double> SlimTreed;
typedef SlimTree<float>  SlimTreef;

#endif // _SLIMTREE_HXX
