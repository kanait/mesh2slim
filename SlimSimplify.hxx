////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIMSIMPLIFY_HXX
#define _SLIMSIMPLIFY_HXX 1

#include "pq.h"
#include "MeshRApp.hxx"
#include "SlimBall.hxx"
#include "SlimTreeApp.hxx"

#include "MeshRDG.hxx"

#define SLIM_ERROR 0
#define SLIM_LEVEL 1

#define CREASE_WEIGHT 10000.0
// #define CREASE_WEIGHT 1000.0

class SlimSimplify : public MeshRdApp, public SlimTreedApp {

public:
  
  SlimSimplify( MeshRd& mesh, SlimTreed& slimtree ) : construct_type_(SLIM_LEVEL),
						      isStorePIEM_(false),
						      isStoreSlimBall_(false),
						      isUpdate_(false) {
    setMesh( mesh );
    setSlimTree( slimtree );
  };
  ~SlimSimplify(){};

  void apply();
  void initialize();
  void createNodePIEMs();
  void setPQ();
  void simplify();
  void polyToSlimBall();
  double optimize( DGEdge*, bool );
  void update( DGEdge* );
  void toSlimTree();

  void clear() {

    pq_.clear();
    dg_.clear();

    slimballs_.clear();

  };

  SlimBalld* createSlimBall() {
    SlimBalld* sb = new SlimBalld;
    slimballs_.push_back( sb );
    return sb;
  };

  void addSlimBall( SlimBalld* sb ) { 
    slimballs_.push_back( sb );
  };

  void setConstructType( unsigned short f ) { construct_type_ = f; };

  void setIsUpdate( bool f ) { isUpdate_ = f; };
  bool isUpdate() const { return isUpdate_; };

  void setIsStoreSlimBall( bool f ) { isStoreSlimBall_ = f; };
  bool isStoreSlimBall() const { return isStoreSlimBall_; };

  void setIsStorePIEM( bool f ) { isStorePIEM_ = f; };
  bool isStorePIEM() const { return isStorePIEM_; };
  

private:

  // SLIM Tree の構築タイプ (0: 誤差順, 1: レベル順)
  unsigned short construct_type_;

  // SLIM balls created
  std::vector<SlimBalld*> slimballs_;

  // Dual Graph
  MeshRDGd dg_;

  // priority queue
  PriorityQueue<PQNoded> pq_;

  // update flag to store slimball
  bool isUpdate_;

  // store PIEM during simplification
  bool isStorePIEM_;

  // store SlimBall during simplification
  bool isStoreSlimBall_;

};

#endif // _SLIMSIMPLIFY_HXX

