////////////////////////////////////////////////////////////////////
//
// $Id: PIEM.hxx 2021/06/05 12:39:35 kanai Exp $
//
// Polygon-Implicit Error Metric class
//
// Copyright (c) 2021 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _PIEM_HXX
#define _PIEM_HXX 1

#include <vector>
using namespace std;

#include "Point3.h"
#include "Vector3.h"
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#define NUM_VEC 11
#define NUM_ELEMENTS 66

class PIEM {

public:
  
  PIEM() : A_(NULL), weightDis_(1.0), weightNrm_(1.0) {
    resize( NUM_ELEMENTS );
    init();
  };
  
  PIEM( PIEM& p ) : A_(NULL), weightDis_(1.0), weightNrm_(1.0) {
    resize( NUM_ELEMENTS );
    set( p );
  };

  ~PIEM(){
    clearMatrixVec();
    clear();
  };

  unsigned int sizeDis() const { return elmDis_.size(); };
  unsigned int sizeNorm() const { return elmNorm_.size(); };
  void resize( int n ) { 
    elmDis_.resize( n ); 
    elmNorm_.resize( n ); 
  };
  void init() { 
    for ( int i = 0; i < NUM_ELEMENTS; ++i ) 
      {
	elmDis_[i] = 0.0; 
	elmNorm_[i] = 0.0; 
      }
  };
  void clear() { elmDis_.clear(); elmNorm_.clear(); };

  void set( const PIEM& p ) {
    for ( unsigned int i = 0; i < NUM_ELEMENTS; ++i ) 
      {
	elmDis_[i]  = p.elmDis(i);
	elmNorm_[i] = p.elmNorm(i);
      }
    area_ = p.area();
  };

  void add( const PIEM& q1, const PIEM& q2 ) {
    for ( unsigned int i = 0; i < NUM_ELEMENTS; ++i )
      {
	elmDis_[i]  = q1.elmDis(i)  + q2.elmDis(i);
	elmNorm_[i] = q1.elmNorm(i) + q2.elmNorm(i);
      }
    area_ = q1.area() + q2.area();
  };

  void add( const PIEM& q ) {
    for ( unsigned int i = 0; i < NUM_ELEMENTS; ++i )
      {
	elmDis_[i]  += q.elmDis(i);
	elmNorm_[i] += q.elmNorm(i);
      }
    area_ += q.area();
  };

  double area() const { return area_; };
  void setArea( double f ) { area_ = f; };

  double elmDis( int i ) const { return elmDis_[i]; };
  double elmNorm( int i ) const { return elmNorm_[i]; };

  void setElementDis( int i, double d ) { elmDis_[i] = d; };
  void addElementDis( int i, double d ) { elmDis_[i] += d; };
  void addElementDis( int i, double d, double w ) { elmDis_[i] += (w*d); };
  void setElementNorm( int i, double d ) { elmNorm_[i] = d; };
  void addElementNorm( int i, double d ) { elmNorm_[i] += d; };
  void addElementNorm( int i, double d, double w ) { elmNorm_[i] += (w*d); };

  void addPIEMElements( Point3d& p0, Point3d& p1, Point3d& p2, Vector3d& nm ) {
    addPIEMElementsDisA( p0, p1, p2 );
    addPIEMElementsNormA( p0, p1, p2 );
    addPIEMElementsNormB( p0, p1, p2, nm );
    addPIEMElementsNormC( nm );
  };
  void addPIEMElementsDisA( Point3d&, Point3d&, Point3d& );
  void addPIEMElementsNormA( Point3d&, Point3d&, Point3d& );
  void addPIEMElementsNormB( Point3d&, Point3d&, Point3d&, Vector3d& );
  void addPIEMElementsNormC( Vector3d& );

  // PIEM から A, b, c を作成
  void copyPIEMToMatrixVec();

  void initMatrixVec() {
    A_ = new double*[ NUM_VEC ];
    for ( int i = 1; i < NUM_VEC; ++i ) A_[i] = new double[ NUM_VEC ];
//     A_.resize( NUM_VEC );
//     for ( int i = 0; i < NUM_VEC; ++i ) A_[i].resize( NUM_VEC );
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

  bool optimizeMKL( std::vector<double>& );
  bool optimize( std::vector<double>& );

  void setWeightDis( double w ) { weightDis_ = w; };
  void setWeightNorm( double w ) { weightNrm_ = w; };

#if 0
  void setMatrixVec( Point3d& p0, Point3d& p1, Point3d& p2, Vector3d& n ) {
    setDisA( p0, p1, p2 );
    setNormA( p0, p1, p2 );
    setNormB( p0, p1, p2, n );
    setNormC( n );
  };

  void setDisA( Point3d&, Point3d&, Point3d& );
  void setNormA( Point3d&, Point3d&, Point3d& );
  void setNormB( Point3d&, Point3d&, Point3d&, Vector3d& );
  void setNormC( Vector3d& );
#endif

private:

  // 66 elements x 2
  std::vector<double> elmDis_;
  std::vector<double> elmNorm_;
  double area_;

  // matrix, vector for optimizing implicit parameters
//   std::vector< std::vector<double> > A_;
  double** A_;
  std::vector<double> b_;
  double c_;

  // weight for adjusting scale of two energies
  double weightDis_;
  double weightNrm_;

};

#endif // _PIEM_HXX

