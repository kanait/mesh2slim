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
    
    // $BLL@Q$N9g7W$r7W;;(B
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
    
    // $B8m:94X?t$NMWAG$N7W;;(B
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


  // PIEM $B$NOB;;$G1"4X?t$N78?t$r7W;;(B
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
    
    // $B%5%]!<%HCf?4!$H>7B$N7W;;(B ($BFs$D$N;R6!$NCf?4$HH>7B$K$h$j7W;;(B)
    calcCenterRadius( centerA, centerB, center, 
		      radiusA, radiusB, radius );

    // $B1"4X?t$N78?t$N7W;;(B (PIEM $B$r;H$C$F(B)
    // $B%0%m!<%P%k:BI8$G7W;;(B
    std::vector<double> coeff_global( QUADRATIC_COEFF );
    
    if ( optimizeBasisFunction( coeff_global ) == false )
      {
	cout << "Warning: optimization failed. " << endl;

	// $B78?t$N7W;;$K<:GT$7$?>l9g!$;R6!$N78?t$rB-$7$F#2$G3d$k(B
	std::vector<double> coeffA_global( QUADRATIC_COEFF );
	std::vector<double> coeffB_global( QUADRATIC_COEFF );
	toGlobalCoord( quadA.coeffs(), centerA, coeffA_global );
	toGlobalCoord( quadB.coeffs(), centerB, coeffB_global );
	averageCoeff( coeffA_global, coeffB_global, coeff_global );
      }
    
    double error = piem().error( coeff_global );

    // $B78?t$r6I=j:BI87O$KJQ49(B
    toLocalCoord( coeff_global, center, coeff );
    
    return error;
  }

  // on-the-fly $B$G1"4X?t$N78?t$r7W;;(B
  double optimize( std::vector<unsigned int> face_list, 
		   const Point3d& centerA, const Point3d& centerB, 
		   Point3d& center, 
		   double radiusA, double radiusB, 
		   double* radius,
		   Quadraticd& quadA, Quadraticd& quadB, 
		   std::vector<double>& coeff ) {
    
    if ( !piem_ ) createPIEM();

    // $B%5%]!<%HCf?4!$H>7B$N7W;;(B ($BFs$D$N;R6!$NCf?4$HH>7B$K$h$j7W;;(B)
    calcCenterRadius( centerA, centerB, center, 
		      radiusA, radiusB, radius );
//     calcCenterRadius( face_list, center, radius );

    // $B%0%m!<%P%k:BI8$G7W;;(B
    std::vector<double> coeff_global( QUADRATIC_COEFF );

    if ( optimizeBasisFunction( face_list, coeff_global ) == false )
      {
	cout << "Warning: optimization failed. " << endl;

	// $B78?t$N7W;;$K<:GT$7$?>l9g!$;R6!$N78?t$rB-$7$F#2$G3d$k(B
	std::vector<double> coeffA_global( QUADRATIC_COEFF );
	std::vector<double> coeffB_global( QUADRATIC_COEFF );
	toGlobalCoord( quadA.coeffs(), centerA, coeffA_global );
	toGlobalCoord( quadB.coeffs(), centerB, coeffB_global );
	averageCoeff( coeffA_global, coeffB_global, coeff_global );
      }
	
    
    double error = piem().error( coeff_global );

    // $B78?t$r6I=j:BI87O$KJQ49(B
    toLocalCoord( coeff_global, center, coeff );
    
    return error;
  };

  // $B%5%]!<%HCf?4!$H>7B$N7W;;!J;R6!$NCf?4!$H>7B$h$j!K(B
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


  // $B%5%]!<%HCf?4!$H>7B$N7W;;!JLL%j%9%H$h$j!K(B
  void calcCenterRadius( std::vector<unsigned int>& face_list, 
			 Point3d& center, double* radius ) {

    std::vector<unsigned int>& indices = mesh().indices();
    std::vector<double>& points = mesh().points();
    std::vector<double>& fnormals = mesh().fnormals();

    // $B;03Q7A$N=E?4$NJ?6Q$h$j7W;;(B
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
    
    // center $B$H;03Q7A$N5wN%$N:GBg(B = radius
    
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

  // $B0l$D$NLL$+$i78?t$r7W;;!J%]%j%4%s$NLL$r(B Slimball $B$KJQ49$9$k$H$-$KMxMQ!K(B
  void createBasisFunction( unsigned int face_id, std::vector<double>& coeff ) {
    
    // p0 $B$7$+;H$o$J$$(B
    Point3d p0, p1, p2;
    mesh().getFacePoints( face_id, p0, p1, p2 );

    Vector3d nrm;
    mesh().fnormal( face_id, nrm );
    
    // $B%]%j%4%s$NJ?LL$r1"4X?tI=8=$K$9$k(B
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

  // PIEM $B$G7W;;(B
  // $B78?t(B: $B%0%m!<%P%k:BI8$G7W;;(B
  bool optimizeBasisFunction( std::vector<double>& coeff ) {
    PIEM& pi = piem();
    return pi.optimize( coeff );
  };

  // on-the-fly $B$G7W;;(B
  // $B78?t(B: $B%0%m!<%P%k:BI8$G7W;;(B
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
