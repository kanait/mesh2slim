////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2005 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIMTREEAPP_HXX
#define _SLIMTREEAPP_HXX 1

#include "SlimTree.hxx"

template <typename T>
class SlimTreeApp {

public:

  SlimTreeApp() : slim_(NULL) {};
  SlimTreeApp( SlimTree<T>& slim ) { setSlimTree( slim ); };
  virtual ~SlimTreeApp(){};

  void clear() { if ( slim_ ) deleteSlimTree(); };
  void setSlimTree( SlimTree<T>& slim ) { slim_ = &slim; };
  SlimTree<T>& slimtree() { return *slim_; };
  
  void deleteSlimTree() { delete slim_; slim_ = NULL; };
  bool empty() const { return ( slim_ != NULL ) ? false : true; };

private:

  SlimTree<T>* slim_;

};

typedef SlimTreeApp<double> SlimTreedApp;
typedef SlimTreeApp<float>  SlimTreefApp;

#endif // _SLIMAPP_HXX
