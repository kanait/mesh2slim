////////////////////////////////////////////////////////////////////
//
// $Id: pq.h 2019/02/09 16:03:15 kanai Exp $
//
// priority queue template class
//
// Copyright (c) 2002 Takashi Kanai
//
// This software is released under the MIT License.
// http://opensource.org/licenses/mit-license.php
//
////////////////////////////////////////////////////////////////////

#ifndef _PQ_H
#define _PQ_H 1

#include <vector>
using namespace std;

//
// class Node_ must include the following member functions:
//   - set( Node_& )
//   - swap( Node_& )
//   - id()  // index
//   - key() // key
//

template <class T>
class PQNode {

public:

  PQNode(){};
  PQNode( int id, T error ) {
    set( id, error );
  };
  ~PQNode() {};

  void set( int id, T error ) {
    setID( id );
    setError( error );
  };
  
  // id
  void setID( int id ) { id_ = id; };
  // needed
  int id() const { return id_; };

  // error
  void setError( T err ) { error_ = err; };
  void setKey( T err ) { setError( err ); };

  // needed
  T key() const { return error_; };

  // needed
  void set( PQNode& p ) {
    setID( p.id() );
    setError( p.key() );
  };

  // needed
  void swap( PQNode& p ) {
    PQNode temp_p;
    temp_p.set( *this );
    set( p );
    p.set( temp_p );
  };

  ////////////////////////////////////////////////////////////////////

  friend bool operator < ( const PQNode& p, const PQNode& q ) {
    return (p.error_ < q.error_);
  };

  friend bool operator > ( const PQNode& p, const PQNode& q ) {
    return (p.error_ > q.error_);
  };

  friend bool operator == ( const PQNode& p, const PQNode& q ) {
    return (p.error_ == q.error_);
  };

  friend bool operator != ( const PQNode& p, const PQNode& q ) {
    return (p.error_ != q.error_);
  };

  friend std::ostream& operator <<( std::ostream& o, const PQNode& q ) {
    return o << "( " << q.error_ << " )\n";
  };
  
private:

  int id_;

  T error_;
  
};

typedef PQNode<unsigned int> PQNodeui;
typedef PQNode<int> PQNodei;
typedef PQNode<float> PQNodef;
typedef PQNode<double> PQNoded;

template<class Node_>
class PriorityQueue {
  
  std::vector<Node_> list_;

  // Node_->id() --> id of std::vector<Node_> map
  std::vector<int> map_;
  
  int last_;
  
public:

  PriorityQueue(){ last_ = 0; };
  ~PriorityQueue(){};

  void init( int size ) {
    list_.resize( size );
    map_.resize( size );
  };

  void clear() {
    list_.clear();
    map_.clear();
    last_ = 0; 
  };
  
  int size() const { return list_.size(); };
  int list_size() const { return last_; };

  bool empty() const { return ( last_ == 0 ); };
  bool exist( int id ) const {
    if ( map_[id] >= last_ ) return false;
    return true;
  };
    
  int insert( Node_& node ) {
    
    assert( last_ < size() );

#if 0    
    if ( last_ >= size() )
      {
	std::cout << "PriorityQueue Warning: size over." << std::endl;
	return 0;
      }
#endif
    
    list_[last_].set( node );
    map_[node.id()] = last_;
    last_++;
    
    int pid = goUp( last_-1 );
    return pid;
  };

  void pop() {
    if ( !last_) return;    // empty

    int id = list_[0].id(); // id of deleted node
    map_[id] = last_;       // list id of deleted node to last_;
    
    list_[0].set( list_[last_-1] );
    map_[list_[0].id()] = 0;
    
    last_--;

    goDown( 0 );
  };

  // id: Node_.id()
  int update( int id, Node_& node ) {
    list_[map_[id]].set( node );
    int pid = goUp( map_[id] );
    return goDown( pid );
  };

  // id: Node_.id()
  Node_& node( int id ) { return list_[map_[id]]; };
  Node_& top() { return list_[0]; };

  // for testing
  int list_id( int id ) const { return map_[id]; };

private:  
  
  // adjust to parent
  int goUp( int id ) {
    int c = id;
    while ( c > 0 ) {
      int p = parent_num( c );
      if ( list_[p].key() > list_[c].key() )
	{
	  swap( p, c );
	  c = p;
	}
      else
	{
	  return c;
	}
    }
    return c;
  };

  // adjust to child
  int goDown( int id ) {

    if ( last_ == 1 ) return id;
    
    int p = id;
    int c, i2l, i2r;
    while ( p <= parent_num(last_-1) )
      {
	// child: i2l, i2r
	i2l = 2*p+1; i2r = i2l+1;
	if ( i2r == last_ ) // pqh->last does not exist.
	  c = i2l;
	else if ( list_[i2l].key() < list_[i2r].key() )
	  c = i2l;  
	else
	  c = i2r;

	if ( list_[p].key() > list_[c].key() )
	  swap( p, c );
	else
	  return p;
	
	p = c;
      }
    return p;
    
  };
  
  int parent_num( int child ) const { return (child - 1) / 2; };

  void swap( int a, int b )
  {
    // store map id
    int id_a = list_[a].id(); // -> b
    int id_b = list_[b].id(); // -> a
    
    assert( map_[id_a] == a );
    assert( map_[id_b] == b );
    assert( map_[id_a] < list_size() );
    assert( map_[id_b] < list_size() );
    
    // swap list
    list_[a].swap( list_[b] );
    // re-map id
    map_[id_a] = b;
    map_[id_b] = a;
  };
  
};

#endif //_PQ_H
