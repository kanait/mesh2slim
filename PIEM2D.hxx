////////////////////////////////////////////////////////////////////
//
// $Id$
//
// 2D Polygon-Implicit Error Metric class
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _PIEM2D_HXX
#define _PIEM2D_HXX 1

// (6 + 1)
#define NUM_VEC 7
#define NUM_DIS_ELEMENTS 21
#define NUM_NRM_ELEMENTS 28

class PIEM2D {

public:
  
  PIEM2D() : A_(NULL), weightDis_(1.0), weightNrm_(1.0), length_(1.0) {
    resize();
    init();
  };
  
  PIEM2D( PIEM2D& p ) : A_(NULL), weightDis_(1.0), weightNrm_(1.0), length_(1.0) {
    resize();
    set( p );
  };

  ~PIEM2D(){
    clearMatrixVec();
    clear();
  };

  unsigned int sizeDis() const { return elmDis_.size(); };
  unsigned int sizeNorm() const { return elmNorm_.size(); };
  
  void resize() { 
    elmDis_.resize(  NUM_DIS_ELEMENTS ); 
    elmNorm_.resize( NUM_NRM_ELEMENTS ); 
  };
  
  void init() { 
    for ( int i = 0; i < NUM_DIS_ELEMENTS; ++i ) 
	elmDis_[i] = 0.0; 

    for ( int i = 0; i < NUM_NRM_ELEMENTS; ++i ) 
	elmNorm_[i] = 0.0; 
  };

  void clear() { elmDis_.clear(); elmNorm_.clear(); };

  void set( const PIEM2D& p ) {
    for ( unsigned int i = 0; i < NUM_DIS_ELEMENTS; ++i ) 
	elmDis_[i]  = p.elmDis(i);

    for ( unsigned int i = 0; i < NUM_NRM_ELEMENTS; ++i ) 
	elmNorm_[i] = p.elmNorm(i);

    length_ = p.length();
  };

  void add( const PIEM2D& q1, const PIEM2D& q2 ) {
    for ( unsigned int i = 0; i < NUM_DIS_ELEMENTS; ++i )
	elmDis_[i]  = q1.elmDis(i)  + q2.elmDis(i);

    for ( unsigned int i = 0; i < NUM_NRM_ELEMENTS; ++i )
	elmNorm_[i] = q1.elmNorm(i) + q2.elmNorm(i);

    length_ = q1.length() + q2.length();
  };

  void add( const PIEM2D& q ) {
    for ( unsigned int i = 0; i < NUM_DIS_ELEMENTS; ++i )
	elmDis_[i]  += q.elmDis(i);

    for ( unsigned int i = 0; i < NUM_NRM_ELEMENTS; ++i )
	elmNorm_[i] += q.elmNorm(i);

    length_ += q.length();
  };

  double length() const { return length_; };
  void setLength( double f ) { length_ = f; };

  double elmDis( int i ) const { return elmDis_[i]; };
  double elmNorm( int i ) const { return elmNorm_[i]; };

  void setElementDis( int i, double d ) { elmDis_[i] = d; };
  void addElementDis( int i, double d ) { elmDis_[i] += d; };
  void addElementDis( int i, double d, double w ) { elmDis_[i] += (w*d); };
  void setElementNorm( int i, double d ) { elmNorm_[i] = d; };
  void addElementNorm( int i, double d ) { elmNorm_[i] += d; };
  void addElementNorm( int i, double d, double w ) { elmNorm_[i] += (w*d); };

  void addPIEMElements( Point2d& p0, Point2d& p1, Vector2d& nm ) {
    addPIEMElementsDisA( p0, p1 );
    addPIEMElementsNormA( p0, p1 );
    addPIEMElementsNormB( p0, p1, nm );
    addPIEMElementsNormC( nm );
  };
  void addPIEMElementsDisA( Point2d&, Point2d& );
  void addPIEMElementsNormA( Point2d&, Point2d& );
  void addPIEMElementsNormB( Point2d&, Point2d&, Vector2d& );
  void addPIEMElementsNormC( Vector2d& );

  // PIEM ‚©‚ç A, b, c ‚ðì¬
  void copyPIEMToMatrixVec();

  void initMatrixVec() {
    A_ = new double*[ NUM_VEC ];
    for ( int i = 1; i < NUM_VEC; ++i ) A_[i] = new double[ NUM_VEC ];
    b_.resize( NUM_VEC );
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

  double error(  std::vector<double>& );

  bool optimize( std::vector<double>& );

  void setWeightDis( double w ) { weightDis_ = w; };
  void setWeightNorm( double w ) { weightNrm_ = w; };

private:

  // 28 elements x 2
  std::vector<double> elmDis_;
  std::vector<double> elmNorm_;
  double length_;

  // matrix, vector for optimizing implicit parameters
  double** A_;
  std::vector<double> b_;
  double c_;

  // weight for adjusting scale of two energies
  double weightDis_;
  double weightNrm_;

};

#endif // _PIEM2D_HXX


