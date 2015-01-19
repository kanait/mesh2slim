////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Polygon-Implicit Error Metric Handling class
//
// Copyright (c) 2005 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#ifndef _PIEM2DHANDLER_HXX
#define _PIEM2DHANDLER_HXX 1

class PIEM2DHandler {

public:

  PIEM2DHandler() : piem_(NULL), isSetLengthWeight_(true), isSetSumLength_(true) {};
  ~PIEM2DHandler(){};

  void setPIEM( PIEM2D& piem ) { piem_ = &piem; };
  PIEM2D& piem() { return *piem_; };

  PIEM2D* createPIEM() {
    PIEM2D* piem = new PIEM2D;
    setPIEM( *piem );
    return piem;
  };

  void deletePIEM() { if ( piem_ ) delete piem_; piem_ = NULL; };

  bool isSetSumLength() const { return isSetSumLength_; };
  void setIsSetSumLength( bool f ) { isSetSumLength_ = f; };

  bool isSetLengthWeight() const { return isSetLengthWeight_; };
  void setIsSetLengthWeight( bool f ) { isSetLengthWeight_ = f; };

  double calcLength( const Point2d& p0, const Point2d& p1 ) {
    return p0.distance( p1 );
  };

  void calcNormal( const Point2d& p0, const Point2d& p1, Vector2d& n ) {

    n.set( p0.y - p1.y, p1.x - p0.x );
    n.normalize();
    
  };

  void setPIEMSumLength( PIEM2D& piem, std::vector<Point2d>& points ) {
    
    // 面積の合計を計算
    double sum_length = 0.0;
    for ( int i = 0; i < points.size()-1; ++i )
      {
	Point2d& p0 = points[i];
	Point2d& p1 = points[i+1];
	double length = calcLength( p0, p1 );
	sum_length += ( length * length );
      }

    piem.setLength( sum_length );
  };

  void setPIEMElements( PIEM2D& piem, std::vector<Point2d>& points ) {
    
    // 誤差関数の要素の計算
    for ( int i = 0; i < points.size()-1; ++i )
      {
	Point2d& p0 = points[i];
	Point2d& p1 = points[i+1];

	Vector2d nrm;
	calcNormal( p0, p1, nrm );
	double length = calcLength( p0, p1 );

	if ( isSetLengthWeight_ == true )
	  {
	    piem.setWeightDis( length );
	    piem.setWeightNorm( length );
	  }

	piem.addPIEMElements( p0, p1, nrm );
	
      }
  };

  // on-the-fly で陰関数の係数を計算
  double optimize( std::vector<Point2d>& points, 
		   Point2d& center, 
		   double* radius,
		   std::vector<double>& coeff ) {
    
    if ( !piem_ ) createPIEM();

    // サポート中心，半径の計算
    calcCenterRadius( points, center, radius );

    // グローバル座標で計算
#if 0
    std::vector<double> coeff_global( NUM_VEC - 1 );
    optimizeBasisFunction( points, coeff_global );
    double error = piem().error( coeff_global );

    // 係数を局所座標系に変換
    toLocalCoord( coeff_global, center, coeff );
#else
    optimizeBasisFunction( points, coeff );
    double error = piem().error( coeff );
#endif

    deletePIEM();

    return error;
  };

  // サポート中心，半径の計算（面リストより）
  void calcCenterRadius( std::vector<Point2d>& points, Point2d& center, double* radius ) {

    // 三角形の重心の平均より計算
    center.set( .0, .0 );
    for ( int i = 0; i < points.size(); ++i )
      {
	center += points[i];
      }
    center.scale( 1.0 / (double) points.size() );
    
    // center と三角形の距離の最大 = radius
    
    *radius = 0.0;
    int c = 0;
    for ( int i = 0; i < points.size(); ++i )
      {
	double dis = center.distance( points[i] );
	if ( i )
	  {
	    if ( dis > *radius ) *radius = dis;
	  }
	else
	  {
	    *radius = dis;
	  }
      }
  };
  
  // on-the-fly で計算
  // 係数: グローバル座標で計算
  bool optimizeBasisFunction( std::vector<Point2d>& points,
			      std::vector<double>& coeff ) {
    PIEM2D& pi = piem();

    if ( isSetSumLength_ == true )
      setPIEMSumLength( pi, points );

    setPIEMElements( pi, points );

    return pi.optimize( coeff );
  };

  // conversion to local coordinate
  void toLocalCoord( std::vector<double>& g,
		     const Point2d& c,
		     std::vector<double>& l )
  {
    l[0] = g[0];
    l[1] = g[1];
    l[2] = g[2];
    l[3] = 2.0*g[0]*c.x + g[2]*c.y + g[3];
    l[4] = 2.0*g[1]*c.y + g[2]*c.x + g[4];
    l[5] = g[0]*c.x*c.x + g[1]*c.y*c.y + g[2]*c.x*c.y + c.x*g[3] + c.y*g[4] + g[5];
  };

  void toGlobalCoord( std::vector<double>& l,
		      const Point2d& c,
		      std::vector<double>& g )
  {
    g[0] = l[0];
    g[1] = l[1];
    g[2] = l[2];
    g[3] = - 2.0*l[0]*c.x - l[2]*c.y + l[3];
    g[4] = - 2.0*l[1]*c.y - l[2]*c.x + l[4];
    g[5] = l[0]*c.x*c.x + l[1]*c.y*c.y + l[2]*c.x*c.y - c.x*l[3] - c.y*l[4] + l[5];
  };

private:

  bool isSetLengthWeight_;
  bool isSetSumLength_;

  PIEM2D* piem_;

};

#endif // _PIEM2DHANDLER_HXX
