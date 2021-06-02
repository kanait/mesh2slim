////////////////////////////////////////////////////////////////////
//
// $Id: MeshRDG.hxx 2021/06/02 22:58:56 kanai Exp $
//
// Dual Graph class for MeshR
//
// Dual Graph クラス
//
//  頂点: 面, エッジ: 面のメイト に対応
//  
//
// Copyright (c) 2005-2006 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _MESHRDG_HXX
#define _MESHRDG_HXX 1

#include <vector>
#include <map>
using namespace std;

#include "MeshR.hxx"
#include "MeshRApp.hxx"

#include "PIEM.hxx"
#include "SlimBall.hxx"

class DGNode;

class DGEdge {

public:

  DGEdge() : sn_(NULL), en_(NULL), slimball_(NULL), piem_(NULL), isDeleted_(false), 
	     isCrease_(false) {};
  ~DGEdge() {
    clear(1);
  };
  
  // type: 
  // 0: エッジのアップデートのとき
  // 1: 簡略化に伴う消去のとき
  void clear( int type ) {
    sn_ = NULL;
    en_ = NULL;
    setIsDeleted( true );
    setIsCrease( false );
    
    switch (type) 
      {
      case 0:
	// PIEM はここで消去
	deletePIEM();
	// ここにある slimball は消去
	deleteSlimBall();
	break;
      case 1: // ここでは何も消去しない (Node に引き継がれているはず)
	setSlimBallPIEM( NULL, NULL );
	break;
      }
  };

  void setID( int i ) { id_ = i; };
  int id() const { return id_; };

  void setIsDeleted( bool f ) { isDeleted_ = f; };
  bool isDeleted() const { return isDeleted_; };

  void setIsCrease( bool f ) { isCrease_ = f; };
  bool isCrease() const { return isCrease_; };

  void setSlimBall( SlimBalld* slimball ) { slimball_ = slimball; };
  SlimBalld* slimball() { return slimball_; };
  void deleteSlimBall() { if ( slimball_ ) delete slimball_; slimball_ = NULL; };

  void setPIEM( PIEM* piem ) { piem_ = piem; };
  PIEM* piem() { return piem_; };
  void deletePIEM() { if ( piem_) delete piem_; piem_ = NULL; };

  void setSlimBallPIEM( SlimBalld* slimball, PIEM* piem ) {
    setSlimBall( slimball );
    setPIEM( piem );
  };

  DGNode* sn() const { return sn_; };
  DGNode* en() const { return en_; };
  DGNode* anotherNode( DGNode* node ) const { return ( sn_ == node ) ? en_ : sn_; };

  void setSnode( DGNode* node ) { sn_ = node; };
  void setEnode( DGNode* node ) { en_ = node; };
  void set( DGNode* sn, DGNode* en ) { sn_ = sn; en_ = en; };

private:

  int id_;

  DGNode* sn_;
  DGNode* en_;

  bool isDeleted_;

  // PIEM
  PIEM* piem_;

  // slim ball created
  SlimBalld* slimball_;

  // crease
  bool isCrease_;
};

class DGNode {

public:

  DGNode() : slimball_(NULL), piem_(NULL), isDeleted_(false) {};

  ~DGNode(){ clear(); };

  // 消去は簡略化の実行のときのみ
  void clear() {
    face_list_.clear();
    edges_.clear();
    setIsDeleted( true );

    // slimball の方は消去しない
    slimball_ = NULL;

    // PIEM はここで消去
    if ( piem_ ) 
      {
	delete piem_;
	piem_ = NULL;
      }
  };

  void add( DGNode* node, DGEdge* orig ) {
    
    // face list
    std::vector<unsigned int>& face_list = node->face_list();
    for ( int i = 0; i < face_list.size(); ++i )
      {
	addFaceList( face_list[i] );
      }

    // edges
    std::vector<DGEdge*>& edges = node->edges();
    for ( int i = 0; i < edges.size(); ++i )
      {
	DGEdge* edge = edges[i];

	if ( edge == NULL ) continue;
	if ( edge->isDeleted() == true ) continue;
	if ( edge == orig ) continue;
	
	// edge の両端のノードのうち, en -> sn とする
	if ( edge->sn() == node ) edge->setSnode( this );
	if ( edge->en() == node ) edge->setEnode( this );

	insertEdge( edge );
      }
    deleteRedundantEdges();
    
  };

  void setSlimBall( SlimBalld* slimball ) { slimball_ = slimball; };
  SlimBalld* slimball() { return slimball_; };
  void deleteSlimBall() { if ( slimball_ ) delete slimball_; slimball_ = NULL; };

  void setPIEM( PIEM* piem ) { piem_ = piem; };
  PIEM* piem() { return piem_; };
  void deletePIEM() { if ( piem_) delete piem_; piem_ = NULL; };

  void setSlimBallPIEM( SlimBalld* slimball, PIEM* piem ) {
    setSlimBall( slimball );
    setPIEM( piem );
  };

  void setIsDeleted( bool f ) { isDeleted_ = f; };
  bool isDeleted() const { return isDeleted_; };

  std::vector<unsigned int>& face_list() { return face_list_; };
  void addFaceList( unsigned int id ) {
    face_list_.push_back( id );
  };
  void deleteFaceList() { face_list_.clear(); };

  std::vector<DGEdge*>& edges() { return edges_; };

  void insertEdge( DGEdge* edge ) { edges_.push_back( edge ); };

  DGEdge* findEdge( DGNode* node ) {
    for ( int i = 0; i < edges_.size(); ++i )
      {
	if ( edges_[i]->isDeleted() ) continue;
	if ( (edges_[i]->sn() == this) && (edges_[i]->en() == node) ||
	     (edges_[i]->sn() == node) && (edges_[i]->en() == this) )
	  return edges_[i];
      }
    return NULL;
  };

  void deleteRedundantEdges() {
//     std::cout << "delete redudant edges ... " << std::endl;
    int dcount = 0;
    for ( int i = 0; i < edges_.size(); ++i )
      {
	DGEdge* edges0 = edges_[i];

	if ( edges0->isDeleted() == true ) continue;

	for ( int j = 0; j < edges_.size(); ++j )
	  {
	    DGEdge* edges1 = edges_[j];
	    
	    if ( edges1->isDeleted() == true ) continue;
	    if ( edges0 == edges1 ) continue;
	    
	    if ( ((edges0->sn() == edges1->sn()) && (edges0->en() == edges1->en())) ||
		 ((edges0->sn() == edges1->en()) && (edges0->en() == edges1->sn())) )
	      {
		// edge1 が crease だったら，edge0 に引き継ぐ
		if ( edges1->isCrease() ) edges0->setIsCrease( true );

		edges1->clear( 0 );
		++( dcount );
	      }
	  }
      }
//     std::cout << "done. num. of deleted edges: " << dcount << std::endl;
  };

private:

  // neighbor edges
  std::vector<DGEdge*> edges_;

  // mesh face list
  std::vector<unsigned int> face_list_;

  bool isDeleted_;

  // slim ball created
  SlimBalld* slimball_;

  // PIEM
  PIEM* piem_;

};

template <typename T>
class MeshRDG : public MeshRApp<T> {

public:

  MeshRDG(){ init(); };
  MeshRDG( MeshR<T>& mesh ) { setMesh(mesh); createDG(); init(); };
  ~MeshRDG(){};
  
  void init() { edge_id_ = 0; };

  void clear() {
    for ( int i = 0; i < nodes_.size(); ++i ) 
      {
	if ( nodes_[i] ) delete nodes_[i];
      }
    nodes_.clear();
    for ( int i = 0; i < edges_.size(); ++i ) 
      {
	if ( edges_[i] ) delete edges_[i];
      }
    edges_.clear();
  };


  DGNode* createNode() {
    DGNode* node = new DGNode;
    nodes_.push_back( node );
    node->setIsDeleted( false );
    return node;
  };

  DGEdge* createEdge( DGNode* sn, DGNode* en ) {
    DGEdge* edge = new DGEdge;
    edges_.push_back( edge );
    edge->setIsDeleted( false );
    edge->set( sn, en );
    edge->setID( edge_id_ ); ++( edge_id_ );
    sn->insertEdge( edge );
    en->insertEdge( edge );
    return edge;
  };

  DGNode* node( int i ) { return nodes_[i]; };
  DGEdge* edge( int i ) { return edges_[i]; };

  std::vector<DGNode*>& nodes() { return nodes_; };
  std::vector<DGEdge*>& edges() { return edges_; };

  //
  // indices_[id]   in MeshR
  // vertices_[id0] in MeshRDG 
  // 
  // id / 3 = id0
  //
  // MeshRDG の要素へは 1/3 して格納する
  //
  void createDG() {
    std::cout << "create Dual Graph ..." << endl;

    std::vector<int>& face_mates = MeshRApp<T>::mesh().face_mates();

    // create nodes
    for ( int i = 0; i < face_mates.size(); i += TRIANGLE )
      {
	DGNode* node = createNode();
	
	// ノードに face の ID を登録
	unsigned int face_id = (int) ((float) i / 3.0f);
	node->addFaceList( face_id );
      }

    // create edges
    for ( int i = 0; i < face_mates.size(); i += TRIANGLE )
      {
	unsigned int id = (int) ((float) i / 3.0f);

	unsigned int m[TRIANGLE];
	MeshRApp<T>::mesh().getFaceMates( id, m );
//  	cout << m[0] << " " << m[1] << " " << m[2] << endl;
// 	m[0] = (int) ((float) face_mates[i] / 3.0f);
// 	m[1] = (int) ((float) face_mates[i+1] / 3.0f);
// 	m[2] = (int) ((float) face_mates[i+2] / 3.0f);
	
	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    if ( m[j] == -1 ) continue;
	    DGNode* node = nodes_[ id ];
	    DGNode* mnode = nodes_[ m[j] ];
	    DGEdge* edge;
	    if ( ( edge = node->findEdge( mnode )) == NULL )
	      {
		edge = createEdge( node, mnode );
#if 0
		if ( mesh().dihedralAngle( id, m[j] ) > CREASE_ANGLE ) // 45°
		  {
		    edge->setIsCrease( true );
		  }
#endif
	      }
	  }
      }
    std::cout << "done." << endl;
  };

  //
  // 縮退した面のノード (DGNode) を削除
  //
  void checkDG() {

    std::cout << "delete degenerate nodes (faces) ... " << std::endl;
    
    int count = 0;
    for ( int i = 0; i < nodes_.size(); ++i )
      {
	DGNode* node = nodes_[i];

	// Mesh の面は一つしか登録されていないはず
	std::vector<unsigned int>& face_list = node->face_list();
	unsigned int face_id = face_list[0];
	
// 	cout << "bb" << endl;

	if ( MeshRApp<T>::mesh().checkFaceNormal( face_id ) ) 
	  continue;

	// 縮退している場合, メイトの面同士をくっつける
// 	cout << "aa" << endl;

	std::vector<DGEdge*> edges = node->edges();
	for ( int j = 0; j < edges.size(); ++j )
	  {
	    if ( edges[j]->isDeleted() ) continue;

	    for ( int k = 0; k < edges.size(); ++k )
	      {
		if ( edges[k]->isDeleted() ) continue;
		if ( edges[j] == edges[k] ) continue;

		DGNode* sn = edges[j]->anotherNode( node );
		DGNode* en = edges[k]->anotherNode( node );

		if ( sn->findEdge( en ) == NULL )
		  createEdge( sn, en );

	      }
	  }

// 	cout << "bb" << endl;
	// エッジの消去
	for ( int j = 0; j < edges.size(); ++j )
	  {
	    if ( edges[j]->isDeleted() ) continue;
 	    edges[j]->clear( -1 );
	  }
	    
// 	cout << "cc" << endl;
	    
	// そのノードは消去
	node->clear();
	++count;

// 	cout << "dd" << endl;
	
      }

    std::cout << "delete degenerate nodes (faces) done. " << count << " nodes deleted. " << std::endl;
  };


private:

  // DG node lists
  std::vector<DGNode*> nodes_;
  
  // DG edge lists
  int edge_id_;
  std::vector<DGEdge*> edges_;
  
};

typedef MeshRDG<double> MeshRDGd;
typedef MeshRDG<float>  MeshRDGf;

#endif // _MESHRDG_HXX


