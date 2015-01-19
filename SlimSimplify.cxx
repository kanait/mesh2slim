////////////////////////////////////////////////////////////////////
//
// $Id$
//
// Copyright (c) 2005-2006 by RIKEN. All rights reserved. 
//
////////////////////////////////////////////////////////////////////

#include "MeshRDG.hxx"

#include "timer.hxx"
#include "SlimBall.hxx"
#include "PIEMHandler.hxx"
#include "SlimSimplify.hxx"

void SlimSimplify::apply()
{
  cout << "slim creation begin. " << endl;

  Timer t;
  double time0 = t.get_seconds();

  initialize();
  
  // メッシュの面を陰関数に変換
  polyToSlimBall();

  setPQ();

  simplify();

  toSlimTree();

  double time1 = t.get_seconds();
  cout << "slim creation end. num of slimballs: " << slimballs_.size() << endl;
  cout << "ellapsed time: " << time1 - time0 << " sec. " << endl;
}

void SlimSimplify::initialize()
{
  cout << "initialize ... " << endl;

  // Face の面法線の作成
  mesh().createFaceNormals( true );

  // Face の面の面積の作成
  mesh().createFaceAreas();

  // Face のメイトの作成
  if ( mesh().face_mates().empty() )
    mesh().createFaceMates();

  // Dual Graph の作成
  dg_.setMesh( mesh() );
  dg_.createDG();

  // 縮退した面のノードを削除する
  dg_.checkDG();

  // Face メイトはもういらないから削除
  mesh().clearFaceMates();

  cout << "initialize done. " << endl;
}

//
// ポリゴンの面を Slimball に変換
//
void SlimSimplify::polyToSlimBall()
{
  cout << "convert original polygons to slimballs ... " << endl;
  std::vector<DGNode*>& nodes = dg_.nodes();
  for ( int i = 0; i < nodes.size(); ++i )
    {
      DGNode* node = nodes[i];
      if ( node->isDeleted() ) continue;

      // SlimBall の作成
      SlimBalld* slimball = createSlimBall();
      
      // DGNode に SlimBall を追加
      node->setSlimBall( slimball );

      PIEMHandler piem_handle( mesh() );

      // サポートの中心と半径を計算
      Point3d center;
      double radius;
      std::vector<unsigned int>& face_list = node->face_list();
      piem_handle.calcCenterRadius( face_list, center, &radius );

      // 二次関数の係数を計算
      Quadraticd* quad = new Quadraticd;
      std::vector<double> coeff_global( QUADRATIC_COEFF );
      // Mesh の面は一つしか登録されていないはず
      unsigned int face_id = face_list[0];
      piem_handle.createBasisFunction( face_id, coeff_global );
      piem_handle.toLocalCoord( coeff_global, center, quad->coeffs() );

#if 0
      Point3d p0, p1, p2;
      mesh().getFacePoints( face_id, p0, p1, p2 );
      Point3d p0sub( p0 - center );
      Point3d p1sub( p1 - center );
      Point3d p2sub( p2 - center );
      cout << "eval " << quad->poly( p0sub ) << " " << quad->poly( p1sub ) 
	   << " " << quad->poly( p2sub ) << endl;
#endif

      // 登録
      slimball->setCenter( center );
      slimball->setSupport( radius );
      slimball->addBasisFunction( *quad );

      // 誤差を登録しておく
      slimball->setUserDefined( 0.0 );
    }
  cout << "convert original polygons to slimballs done. " << endl;
}

double SlimSimplify::optimize( DGEdge* edge, bool isUpdate )
{
  if ( edge->slimball() ) edge->deleteSlimBall();
  if ( edge->piem() ) edge->deletePIEM();
  
  // SlimBall の作成
  SlimBalld* slimball = new SlimBalld;
      
  // PIEM の作成
  PIEM* piem = new PIEM;

  // 両端のノードから face list を取得
  std::vector<unsigned int> face_list;
  DGNode* sn = edge->sn();
  std::vector<unsigned int>& sfl = sn->face_list();
  for ( int i = 0; i < sfl.size(); ++i ) { face_list.push_back( sfl[i] ); };
  DGNode* en = edge->en();
  std::vector<unsigned int>& efl = en->face_list();
  for ( int i = 0; i < efl.size(); ++i ) { face_list.push_back( efl[i] ); };

//       cout << "face_list " << face_list.size() << endl;

  // PIEM による最適関数の計算

  PIEMHandler piem_handle( mesh() );
  piem_handle.setPIEM( *piem );
  
  Point3d center;
  double radius;
  Quadraticd* quad = new Quadraticd;

  double error;
  if ( isStorePIEM() == false )
    {
      error = piem_handle.optimize( face_list, 
				    sn->slimball()->center(),
				    en->slimball()->center(),
				    center, 
				    sn->slimball()->support(),
				    en->slimball()->support(),
				    &radius, 
				    (Quadraticd&) sn->slimball()->bf(),
				    (Quadraticd&) en->slimball()->bf(),
				    quad->coeffs() );
    }
  else
    {
      error = piem_handle.optimize( sn->piem(), en->piem(),
				    sn->slimball()->center(),
				    en->slimball()->center(),
				    center, 
				    sn->slimball()->support(),
				    en->slimball()->support(),
				    &radius, 
				    (Quadraticd&) sn->slimball()->bf(),
				    (Quadraticd&) en->slimball()->bf(),
				    quad->coeffs() );
    }

  // crease のときは誤差関数値を 1000 倍にする
//   if ( edge->isCrease() == true ) error *= 1000.0;
  if ( edge->isCrease() == true ) error *= CREASE_WEIGHT;
      
  slimball->setCenter( center );
  slimball->setSupport( radius );
  slimball->addBasisFunction( *quad );
  // 誤差を登録しておく
  slimball->setUserDefined( error );

  // update の中で使われる時以外は slimball を消去する
  if ( isUpdate == false )
    {
      delete slimball; slimball = NULL;
    }

//   if ( (isStorePIEM() == false) && (isUpdate == false) )
  if ( isUpdate == false )
    {
      delete piem; piem = NULL;
    }

  edge->setSlimBallPIEM( slimball, piem );
  
  return error;
}

// 頂点のノードに PIEM を作成
void SlimSimplify::createNodePIEMs()
{
  cout << "create node PIEMs ... " << endl;
  PIEMHandler piem_handle( mesh() );
  
  for ( int i = 0; i < dg_.nodes().size(); ++i )
    {
      DGNode* node = dg_.node(i);
      if ( node->isDeleted() == true ) continue;
	    
      SlimBalld* slimball = node->slimball();
      node->setPIEM( piem_handle.createPIEM( node->face_list() ) );

      // ノードの face_list はもういらない（はず）なので消去
      node->deleteFaceList();
    }
  cout << "create node PIEMs done. " << endl;
}

void SlimSimplify::setPQ()
{
  cout << "create priority queue ... " << endl;

  // queue の初期化
  pq_.init( dg_.edges().size() );

//   cout << " isStorePIEM " << isStorePIEM() << endl;
  if ( isStorePIEM() )
    {
      createNodePIEMs();
    }

  // MeshDG のエッジから PQ のノードを作成
  for ( int k = 0; k < dg_.edges().size(); ++k )
    {
      DGEdge* edge = dg_.edge(k);
      if ( edge->isDeleted() ) continue;
      
      double error = optimize( edge, false );
      
//       cout << "edge " << k << " error " << error << endl;

      // PQ ノードの作成
      PQNoded pqnode( k, error );
      pq_.insert( pqnode );
    }

  cout << "create priority queue done. " << endl;
}

void SlimSimplify::update( DGEdge* edge )
{
//   if ( isStoreSlimBall() == false )
//     {
  optimize( edge, true );
//     }

  // slimball をリストに追加
  addSlimBall( edge->slimball() );

//   edge->slimball()->print();

  // SlimBall の親子関係の構築
  DGNode* sn = edge->sn();
  if ( sn->slimball() )
    {
      sn->slimball()->setParent( edge->slimball() );
      edge->slimball()->addChild( sn->slimball() );
    }

  DGNode* en = edge->en();
  if ( en->slimball() )
    {
      en->slimball()->setParent( edge->slimball() );
      edge->slimball()->addChild( en->slimball() );
    }

  // エッジの slimball と PIEM を sn に引き継ぐ
  // sn にもともと付いている PIEM は消去, slimball はリンクを切る
  
  sn->setSlimBall( edge->slimball() );
  
  sn->deletePIEM();
  sn->setPIEM( edge->piem() );

  // 引継ぎ終わったので NULL にしておく
  edge->setSlimBall( NULL );
  edge->setPIEM( NULL );

  // en の face_list と edges を sn に追加
  sn->add( en, edge );

//   cout << "\t new node ... num of faces: " << sn->face_list().size() << endl;

  // en と edge を clear
  en->clear();
  edge->clear( 1 );
}

void SlimSimplify::simplify()
{
  cout << "simplification ... " << endl;

  bool changeNodeStoreType = false;
  int num_faces = mesh().numFaces();

  int count = 0;
  while ( !(pq_.empty()) )
    {
      PQNoded& node = pq_.top();
      int min_id = node.id();
      double min_error = node.key();
      pq_.pop();

      DGEdge* edge = dg_.edge( min_id );
      
      if ( edge->isDeleted() == true ) continue;

//       if ( edge->isCrease() )
// 	cout << "edge isCrease" << edge->isCrease() << endl;
      
      ++count;
      if ( !(count % 1000) )
	cout << "count " << count << "/" << num_faces << " min_id " << min_id << " min_error " << min_error << endl;

      // 処理すべきノードのサイズが もとの面のサイズの 1/10 になったら
      // PIEM 保存に切り替える
//       cout << "tmp " << slimballs_.size() - mesh().numFaces() << endl;
      if ( ( count > num_faces * .9 ) && ( changeNodeStoreType == false ) && ( isStorePIEM() == false) )
	{
	  createNodePIEMs();
	  setIsStorePIEM( true );
	  changeNodeStoreType = true;
	}

      DGNode* sn = edge->sn();

      update( edge );

      // ノードに隣接するエッジの slimball と PIEM をアップデート
      std::vector<DGEdge*> edges = sn->edges();
      for ( int i = 0; i < edges.size(); ++i )
	{
	  DGEdge* ne = edges[i];
	  if ( ne->isDeleted() ) continue;

	  // slimball と PIEM をアップデート
	  double error = optimize( ne, false );
	  PQNoded pqnode( ne->id(), error );
	  pq_.update( ne->id(), pqnode );
	}
    }

  cout << "simplification done. " << endl;
}

// slimball の並び替え
void SlimSimplify::toSlimTree()
{
  cout << "construct SlimTree ... " << endl;
  if ( construct_type_ == SLIM_ERROR )
    {
      int c = 0; 
      for ( int i = slimballs_.size() - 1; i >= 0; --i )
	{
	  slimtree().setSlimBall( slimballs_[i] ); ++c;
	}
    }
  else if ( construct_type_ == SLIM_LEVEL )
    {
      // level 0 の slimball を取得
      //std::vector<SlimBall*> slim_level;

      int level = 0;
      int c = 0;
      for ( int i = slimballs_.size() - 1; i >= 0; --i )
	{
	  if ( !( slimballs_[i]->parent() ) )
	    {
	      slimtree().setSlimBall( slimballs_[i] ); 

	      ++c;
	    }
	}
      
      ++( level );

      int start = 0;
      int end = c;
      while ( 1 ) 
	{
	  for ( int i = start; i < end; ++i )
	    {
	      SlimBalld* slimball = slimtree().slimball(i);
	      if ( !slimball ) continue;

	      std::vector<SlimBalld*>& childs = slimball->childs();
	      for ( int j = 0; j < childs.size(); ++j )
		{
		  // 葉ノードを保存しない場合は次の行をコメントアウトする
// 		  if ( childs[j]->childs().size() )

		  slimtree().setSlimBall( childs[j] ); 

		  ++c;
		}
	    }
	  
	  start = end;
	  end = c;
	  
	  if ( end == slimballs_.size() ) break;

	  ++( level );
	}
	  
    }
  cout << "construct SlimTree done. " << endl;
}
