////////////////////////////////////////////////////////////////////
//
// $Id$
//
// 2D Polygon-Implicit Error Metric class
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _POINTFIT2D_HXX
#define _POINTFIT2D_HXX 1

// (6 + 1)
#define NUM_VEC 7

class PointFit2D {

public:

  PointFit2D() : A_(NULL), weightDis_(1.0), weightNrm_(1.0) {
    initMatrixVec();
  };
  
  ~PointFit2D() {
    clearMatrixVec();
  };

  void setWeightDis( double w ) { weightDis_ = w; };
  void setWeightNorm( double w ) { weightNrm_ = w; };

  double error(  std::vector<double>& );

  bool optimize( std::vector<double>& );

  void initMatrixVec() {
    A_ = new double*[ NUM_VEC ];
    for ( int i = 1; i < NUM_VEC; ++i ) 
      {
	A_[i] = new double[ NUM_VEC ];
	for ( int j = 1; j < NUM_VEC; ++j )
	  {
	    A_[i][j] = 0.0;
	  }
      }
    
    b_.resize( NUM_VEC );
    for ( int i = 1; i < NUM_VEC; ++i ) 
      {
	b_[i] = 0.0;
      }

    c_ = 0.0;
    
  };

  void clearMatrixVec() {
    
    if ( A_ )
      {
	for ( int i = 1; i < NUM_VEC; ++i ) delete A_[i];
	delete A_;
	A_ = NULL;
      }

    b_.clear();
  };

  void addElementDis( int i, int j, double d ) { A_[i+1][j+1] += (weightDis_*d); };
  void addElementNormA( int i, int j, double d ) { A_[i+1][j+1] += (weightNrm_*d); };
  void addElementNormb( int i, double d ) { b_[i+1] += (weightNrm_*d); };
  void addElementNormc( double d ) { c_ += (weightNrm_*d); };

  void addPointFitElements( Point2d& p, Vector2d& nm ) {
    addPointFitElementsDisA( p );
    addPointFitElementsNormA( p );
    addPointFitElementsNormB( p, nm );
    addPointFitElementsNormC( nm );
  };
  void addPointFitElementsDisA( Point2d& );
  void addPointFitElementsNormA( Point2d& );
  void addPointFitElementsNormB( Point2d&, Vector2d& );
  void addPointFitElementsNormC( Vector2d& );

private:

  // matrix, vector for optimizing implicit parameters
  double** A_;
  std::vector<double> b_;
  double c_;

  // weight for adjusting scale of two energies
  double weightDis_;
  double weightNrm_;

};

#endif // _POINTFIT2D_HXX

