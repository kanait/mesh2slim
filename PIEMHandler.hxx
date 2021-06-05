////////////////////////////////////////////////////////////////////
//
// $Id: PIEMHandler.hxx 2021/06/05 12:40:00 kanai Exp $
//
// Polygon-Implicit Error Metric Handling class
//
// Copyright (c) 2021 Takashi Kanai
// Released under the MIT license
//
////////////////////////////////////////////////////////////////////

#ifndef _PIEMHANDLER_HXX
#define _PIEMHANDLER_HXX 1

#include <vector>
using namespace std;

#include "Point3.h"
#include "Vector3.h"
#ifdef VM_INCLUDE_NAMESPACE
using namespace kh_vecmath;
#endif // VM_INCLUDE_NAMESPACE

#include "MeshR.hxx"
#include "MeshRApp.hxx"
#include "PIEM.hxx"

class PIEMHandler : public MeshRdApp {

public:

  PIEMHandler() : piem_(NULL) {};
  PIEMHandler( MeshRd& mesh ) : piem_(NULL) { setMesh( mesh ); };
  PIEMHandler( MeshRd& mesh, PIEM& piem ) { setMesh( mesh ); setPIEM( piem ); };
  
  ~PIEMHandler(){};

  void setPIEM( PIEM& piem ) { piem_ = &piem; };
  PIEM& piem() { return *piem_; };

  PIEM* createPIEM() {
    PIEM* piem = new PIEM;
    setPIEM( *piem );
    return piem;
  };

  double error( std::vector<double>& );

  void setPIEMSumArea( PIEM& piem, std::vector<unsigned int> face_list ) {
    
    // 面積の合計を計算
    double sum_area = 0.0;
    for ( int i = 0; i < face_list.size(); ++i )
      {
	unsigned int face_id = face_list[i];
	// face area as a weight
	sum_area += mesh().farea( face_id );
      }

    piem.setArea( sum_area );
  };

  void setPIEMElements( PIEM& piem, std::vector<unsigned int>& face_list ) {
    
    // 誤差関数の要素の計算
    for ( int i = 0; i < face_list.size(); ++i )
      {
	unsigned int face_id = face_list[i];
	
	Vector3d nrm;
	mesh().fnormal( face_id, nrm );
	double area = mesh().farea( face_id );

	Point3d p0, p1, p2;
	mesh().getFacePoints( face_id, p0, p1, p2 );

 	piem.setWeightDis( area );
   	piem.setWeightNorm( area );
	piem.addPIEMElements( p0, p1, p2, nrm );

      }
  };

  PIEM* createPIEM( std::vector<unsigned int>& face_list ) {
    
    PIEM* piem = new PIEM;

    setPIEMSumArea( *piem, face_list );
    setPIEMElements( *piem, face_list );

    return piem;
  };


  // PIEM の和算で陰関数の係数を計算
  double optimize( PIEM* piemA, PIEM* piemB, 
		   const Point3d& centerA, const Point3d& centerB, 
		   Point3d& center, 
		   double radiusA, double radiusB, 
		   double* radius,
		   Quadraticd& quadA, Quadraticd& quadB, 
		   std::vector<double>& coeff )
  {
    if ( !piem_ ) createPIEM();

    piem().add( *piemA, *piemB );
    
    // サポート中心，半径の計算 (二つの子供の中心と半径により計算)
    calcCenterRadius( centerA, centerB, center, 
		      radiusA, radiusB, radius );

    // 陰関数の係数の計算 (PIEM を使って)
    // グローバル座標で計算
    std::vector<double> coeff_global( QUADRATIC_COEFF );
    
    if ( optimizeBasisFunction( coeff_global ) == false )
      {
	cout << "Warning: optimization failed. " << endl;

	// 係数の計算に失敗した場合，子供の係数を足して２で割る
	std::vector<double> coeffA_global( QUADRATIC_COEFF );
	std::vector<double> coeffB_global( QUADRATIC_COEFF );
	toGlobalCoord( quadA.coeffs(), centerA, coeffA_global );
	toGlobalCoord( quadB.coeffs(), centerB, coeffB_global );
	averageCoeff( coeffA_global, coeffB_global, coeff_global );
      }
    
    double error = piem().error( coeff_global );

    // 係数を局所座標系に変換
    toLocalCoord( coeff_global, center, coeff );
    
    return error;
  }

  // on-the-fly で陰関数の係数を計算
  double optimize( std::vector<unsigned int> face_list, 
		   const Point3d& centerA, const Point3d& centerB, 
		   Point3d& center, 
		   double radiusA, double radiusB, 
		   double* radius,
		   Quadraticd& quadA, Quadraticd& quadB, 
		   std::vector<double>& coeff ) {
    
    if ( !piem_ ) createPIEM();

    // サポート中心，半径の計算 (二つの子供の中心と半径により計算)
    calcCenterRadius( centerA, centerB, center, 
		      radiusA, radiusB, radius );
//     calcCenterRadius( face_list, center, radius );

    // グローバル座標で計算
    std::vector<double> coeff_global( QUADRATIC_COEFF );

    if ( optimizeBasisFunction( face_list, coeff_global ) == false )
      {
	cout << "Warning: optimization failed. " << endl;

	// 係数の計算に失敗した場合，子供の係数を足して２で割る
	std::vector<double> coeffA_global( QUADRATIC_COEFF );
	std::vector<double> coeffB_global( QUADRATIC_COEFF );
	toGlobalCoord( quadA.coeffs(), centerA, coeffA_global );
	toGlobalCoord( quadB.coeffs(), centerB, coeffB_global );
	averageCoeff( coeffA_global, coeffB_global, coeff_global );
      }
	
    
    double error = piem().error( coeff_global );

    // 係数を局所座標系に変換
    toLocalCoord( coeff_global, center, coeff );
    
    return error;
  };

  // サポート中心，半径の計算（子供の中心，半径より）
  void calcCenterRadius( const Point3d& centerA, const Point3d& centerB, 
			 Point3d& center, 
			 double radiusA, double radiusB, 
			 double* radius )
  {
//     center.set( centerA + centerB );
//     center.scale( .5 );

//     double rad = ( radiusA > radiusB ) ? radiusA : radiusB;
//     *radius = rad + centerA.distance( centerB ) * .5;

    double dis = centerA.distance( centerB );

    if( (radiusA > radiusB) && ((dis+radiusB) < radiusA) )
      {
	*radius = radiusA;
	center.set( centerA );
      }
    else if( (radiusB > radiusA) && ((dis+radiusA) < radiusB) )
      {
	*radius = radiusB;
	center.set( centerB );
      }
    else
      {
        *radius = 0.5 * (radiusA+radiusB+dis);
        double w1 = (dis+radiusA-radiusB)/(2.0*dis);
        double w2 = (dis-radiusA+radiusB)/(2.0*dis);
	center.set( w1*centerA.x + w2*centerB.x,
		    w1*centerA.y + w2*centerB.y,
		    w1*centerA.z + w2*centerB.z );
      }
  };


  // サポート中心，半径の計算（面リストより）
  void calcCenterRadius( std::vector<unsigned int>& face_list, 
			 Point3d& center, double* radius ) {

    std::vector<unsigned int>& indices = mesh().indices();
    std::vector<double>& points = mesh().points();
    std::vector<double>& fnormals = mesh().fnormals();

    // 三角形の重心の平均より計算
    center.set( .0, .0, .0 );
    for ( int i = 0; i < face_list.size(); ++i )
      {
	unsigned int face_id = face_list[i];
	
	Point3d p0, p1, p2;
	mesh().getFacePoints( face_id, p0, p1, p2 );

	Point3d barycenter( p0 + p1 + p2 ); barycenter.scale( 1.0 / (double) TRIANGLE );
	center += barycenter;
      }
    center.scale( 1.0 / (double) face_list.size() );
    
    // center と三角形の距離の最大 = radius
    
    *radius = 0.0;
    int c = 0;
    for ( int i = 0; i < face_list.size(); ++i )
      {
	unsigned int face_id = face_list[i];
	for ( int j = 0; j < TRIANGLE; ++j )
	  {
	    Point3d p( points[ nXYZ * indices[ TRIANGLE * face_id + j ] ],
		       points[ nXYZ * indices[ TRIANGLE * face_id + j ] + 1 ],
		       points[ nXYZ * indices[ TRIANGLE * face_id + j ] + 2 ] );
	    double dis = center.distance( p );
	    if ( c )
	      {
		if ( dis > *radius ) *radius = dis;
	      }
	    else
	      {
		*radius = dis;
	      }
	    ++c;
	  }
      }
  };
  
  void calcCenterRadius( Point3d& center, double* radius ) {

    std::vector<unsigned int>& indices = mesh().indices();
    std::vector<double>& points = mesh().points();
    std::vector<double>& fnormals = mesh().fnormals();

    // calculate support center and radius

    center.set( .0, .0, .0 );
    for ( int i = 0; i < mesh().numPoints(); ++i )
      {
	center.x += points[ nXYZ * i ];
	center.y += points[ nXYZ * i + 1 ];
	center.z += points[ nXYZ * i + 2 ];
      }
    center.scale( 1.0 / (double) mesh().numPoints() );
//     cout << "center = " << center << endl;
  
    *radius = 0.0;
    for ( int i = 0; i < mesh().numPoints(); ++i )
      {
	double dis = center.distance( Point3d( points[ nXYZ * i ],
					       points[ nXYZ * i + 1 ],
					       points[ nXYZ * i + 2 ] ) );
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

  // 一つの面から係数を計算（ポリゴンの面を Slimball に変換するときに利用）
  void createBasisFunction( unsigned int face_id, std::vector<double>& coeff ) {
    
    // p0 しか使わない
    Point3d p0, p1, p2;
    mesh().getFacePoints( face_id, p0, p1, p2 );

    Vector3d nrm;
    mesh().fnormal( face_id, nrm );
    
    // ポリゴンの平面を陰関数表現にする
    coeff[0] = 0.0;
    coeff[1] = 0.0;
    coeff[2] = 0.0;
    coeff[3] = 0.0;
    coeff[4] = 0.0;
    coeff[5] = 0.0;
    coeff[6] = nrm.x;
    coeff[7] = nrm.y;
    coeff[8] = nrm.z;
    coeff[9] = -( nrm.x * p0.x + nrm.y * p0.y + nrm.z * p0.z );

  };

  // PIEM で計算
  // 係数: グローバル座標で計算
  bool optimizeBasisFunction( std::vector<double>& coeff ) {
    PIEM& pi = piem();
    return pi.optimize( coeff );
  };

  // on-the-fly で計算
  // 係数: グローバル座標で計算
  bool optimizeBasisFunction( std::vector<unsigned int>& face_list,
			      std::vector<double>& coeff ) {
    PIEM& pi = piem();

    setPIEMSumArea( pi, face_list );
    setPIEMElements( pi, face_list );

//     return pi.optimizeMKL( coeff );
    return pi.optimize( coeff );
  };

  void optimizeBasisFunction( const Point3d& center, std::vector<double>& coeff ) {

    std::vector<unsigned int>& indices = mesh().indices();
    std::vector<double>& points = mesh().points();
    std::vector<double>& fnormals = mesh().fnormals();

    if ( !piem_ ) createPIEM();
    PIEM& pi = piem();

    for ( int i = 0; i < mesh().numFaces(); ++i )
      {
	Point3d p0( points[ nXYZ * indices[ TRIANGLE * i ] ],
		    points[ nXYZ * indices[ TRIANGLE * i ] + 1 ],
		    points[ nXYZ * indices[ TRIANGLE * i ] + 2 ] );
	Point3d p1( points[ nXYZ * indices[ TRIANGLE * i + 1 ] ],
		    points[ nXYZ * indices[ TRIANGLE * i + 1 ] + 1 ],
		    points[ nXYZ * indices[ TRIANGLE * i + 1 ] + 2 ] );
	Point3d p2( points[ nXYZ * indices[ TRIANGLE * i + 2 ] ],
		    points[ nXYZ * indices[ TRIANGLE * i + 2 ] + 1 ],
		    points[ nXYZ * indices[ TRIANGLE * i + 2 ] + 2 ] );

	// convert local coordinates
	p0 -= center;
	p1 -= center;
	p2 -= center;

#if 1
	// normal vector (not normalized)
	Vector3d sub1( p1 - p0 );
	Vector3d sub2( p2 - p0 );
	Vector3d nrm;
	nrm.cross( sub1, sub2 );
#endif

      // from face normal
#if 0
	Vector3d nrm( fnormals[ nXYZ * i ],
		      fnormals[ nXYZ * i + 1 ],
		      fnormals[ nXYZ * i + 2 ] );
#endif

	// face area as a weight
	double area = mesh().calcFaceArea( p0, p1, p2 );

	pi.setWeightDis( area );
	pi.setWeightNorm( area );
	pi.addPIEMElements( p0, p1, p2, nrm );

      }

//     pi.optimizeMKL( coeff );
    pi.optimize( coeff );
  };

  // conversion to local coordinate
  void toLocalCoord( std::vector<double>& g,
		     const Point3d& c,
		     std::vector<double>& l )
  {
    l[0] = g[0];
    l[1] = g[1];
    l[2] = g[2];
    l[3] = g[3];
    l[4] = g[4];
    l[5] = g[5];
    l[6] = 2.0 * g[0] * c.x + c.y * g[3] + c.z * g[5] + g[6];
    l[7] = 2.0 * g[1] * c.y + c.x * g[3] + c.z * g[4] + g[7];
    l[8] = 2.0 * g[2] * c.z + c.y * g[4] + c.x * g[5] + g[8];
    l[9] = ( g[0] * c.x * c.x + g[1] * c.y * c.y + g[2] * c.z * c.z +
	     c.x * c.y * g[3] + c.y * c.z * g[4] + c.x * c.z * g[5] + 
	     c.x * g[6] + c.y * g[7] + c.z * g[8] + g[9] );
  };

  void toGlobalCoord( std::vector<double>& l,
		      const Point3d& c,
		      std::vector<double>& g )
  {
    g[0] = l[0];
    g[1] = l[1];
    g[2] = l[2];
    g[3] = l[3];
    g[4] = l[4];
    g[5] = l[5];
    g[6] = ( -2.0 * l[0] * c.x - c.y * l[3] - c.z * l[5] + l[6] );
    g[7] = ( -2.0 * l[1] * c.y - c.x * l[3] - c.z * l[4] + l[7] );
    g[8] = ( -2.0 * l[2] * c.z - c.y * l[4] - c.x * l[5] + l[8] );
    g[9] = ( l[0] * c.x * c.x + l[1] * c.y * c.y + l[2] * c.z * c.z + 
	     c.x * c.y * l[3] + c.y * c.z * l[4] + c.x * c.z * l[5] - 
	     c.x * l[6] - c.y * l[7] - c.z * l[8] + l[9] );
  };

  void averageCoeff( std::vector<double>& cA, std::vector<double>& cB, std::vector<double>& c )
  {
    for ( int i = 0; i < QUADRATIC_COEFF; ++i )
      c[i] = ( cA[i] + cB[i] ) * .5;
  };

private:

  PIEM* piem_;

};

#endif // _PIEMHANDLER_HXX
