////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2005 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _MESHRAPP_HXX
#define _MESHRAPP_HXX 1

#include "MeshR.hxx"

template <class T>
class MeshRApp {

public:

  MeshRApp() : mesh_(NULL) {};
  MeshRApp( MeshR<T>& mesh ) { setMesh( mesh ); };
  virtual ~MeshRApp(){};

  void clear() { if ( mesh_ ) deleteMesh(); };
  void setMesh( MeshR<T>& mesh ) { mesh_ = &mesh; };
  MeshR<T>& mesh() { return *mesh_; };
  
  void deleteMesh() { delete mesh_; mesh_ = NULL; };
  bool empty() const { return ( mesh_ != NULL ) ? false : true; };

private:

  MeshR<T>* mesh_;

};

typedef MeshRApp<double> MeshRdApp;
typedef MeshRApp<float>  MeshRfApp;

#endif // _MESHRAPP_HXX


