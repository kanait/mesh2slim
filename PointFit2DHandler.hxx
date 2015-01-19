////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Polygon-Implicit Error Metric Handling class
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _POINTFIT2DHANDLER_HXX
#define _POINTFIT2DHANDLER_HXX 1

class PointFit2DHandler {

public:

  PointFit2DHandler() : pf_(NULL), isSetLengthWeight_(true), isSetSumLength_(true) {};
  ~PointFit2DHandler(){};

  void setPointFit( PointFit2D& pf ) { pf_ = &pf; };
  PointFit2D& pointfit() { return *pf_; };

  PointFit2D* createPointFit() {
    PointFit2D* pf = new PointFit2D;
    setPointFit( *pf );
    return pf;
  };

  void deletePointFit() { if ( pf_ ) delete pf_; pf_ = NULL; };

  bool isSetSumLength() const { return isSetSumLength_; };
  void setIsSetSumLength( bool f ) { isSetSumLength_ = f; };

  bool isSetLengthWeight() const { return isSetLengthWeight_; };
  void setIsSetLengthWeight( bool f ) { isSetLengthWeight_ = f; };

  void setLength( double length ) { length_ = length; };
  double length() const { return length_; };
  
  void setPointFitElements( PointFit2D& pf, std::vector<Point2d>& points, 
			    std::vector<Vector2d>& normals ) {
    
    pf.setWeightNorm( length() );
	
    // 誤差関数の要素の計算
    for ( int i = 0; i < points.size(); ++i )
      {
	pf.addPointFitElements( points[i], normals[i] );
      }
  };

  // on-the-fly で陰関数の係数を計算
  double optimize( std::vector<Point2d>& points, 
		   std::vector<Vector2d>& normals, 
		   double length,
		   std::vector<double>& coeff ) {
    
    if ( !pf_ ) createPointFit();

    setLength( length );

    optimizeBasisFunction( points, normals, coeff );
    double error = pointfit().error( coeff );

    deletePointFit();

    return error;
  };

  // on-the-fly で計算
  // 係数: グローバル座標で計算
  bool optimizeBasisFunction( std::vector<Point2d>& points,
			      std::vector<Vector2d>& normals, 
			      std::vector<double>& coeff ) {
    PointFit2D& pf = pointfit();

    setPointFitElements( pf, points, normals );

    return pf.optimize( coeff );
  };

private:

  bool isSetLengthWeight_;
  bool isSetSumLength_;

  double length_;

  PointFit2D* pf_;

};

#endif // _POINTFIT2DHANDLER_HXX
