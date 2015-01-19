////////////////////////////////////////////////////////////////////
//
// $Id: $
//
// Copyright (c) 2005 by Takashi Kanai. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _SLIMAPP_HXX
#define _SLIMAPP_HXX 1

#include "Slim.hxx"

template <typename T>
class SlimApp {

public:

  SlimApp() : slim_(NULL) {};
  SlimApp( Slim<T>& slim ) { setSlim( slim ); };
  virtual ~SlimApp(){};

  void clear() { if ( slim_ ) deleteSlim(); };
  void setSlim( Slim<T>& slim ) { slim_ = &slim; };
  Slim<T>& slim() { return *slim_; };
  
  void deleteSlim() { delete slim_; slim_ = NULL; };
  bool empty() const { return ( slim_ != NULL ) ? false : true; };

private:

  Slim<T>* slim_;

};

typedef SlimApp<double> SlimdApp;
typedef SlimApp<float>  SlimfApp;

#endif // _SLIMAPP_HXX
