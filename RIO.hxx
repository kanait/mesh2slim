////////////////////////////////////////////////////////////////////
//
// $Id: RIO.hxx 2021/06/05 12:40:20 kanai Exp $
//
// Copyright (c) 2021 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _RDIO_HXX
#define _RDIO_HXX 1

#include "MeshR.hxx"

template <typename T>
class RIO {

  MeshR<T>* mesh_;

  bool saveTexcoord_;
  bool saveColor_;
  bool saveNormal_;
  bool saveBLoop_;

public:

  RIO() { init(); };
  RIO( MeshR<T>& mesh ) { init(); setMesh( mesh ); };
  virtual ~RIO() {};
  
  void init() {
    mesh_ = NULL;
    saveTexcoord_ = false;
    saveColor_ = false;
    saveNormal_ = false;
    saveBLoop_ = false;
  };
  
  virtual MeshR<T>& mesh() { return *mesh_; };

  bool empty() const { return ( mesh_ == NULL ) ? true : false; };
  void setMesh( MeshR<T>& mesh ) { mesh_ = &mesh; };

  bool isSaveTexcoord() const { return saveTexcoord_; };
  void setSaveTexcoord( bool f ) { saveTexcoord_ = f; };

  bool isSaveColor() const { return saveColor_; };
  void setSaveColor( bool f ) { saveColor_ = f; };

  bool isSaveNormal() const { return saveNormal_; };
  void setSaveNormal( bool f ) { saveNormal_ = f; };

  bool isSaveBLoop() const { return saveBLoop_; };
  void setSaveBLoop( bool f ) { saveBLoop_ = f; };
};

#endif // _RDIO_HXX
