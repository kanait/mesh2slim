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
  
  // ���b�V���̖ʂ��A�֐��ɕϊ�
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

  // Face �̖ʖ@���̍쐬
  mesh().createFaceNormals( true );

  // Face �̖ʂ̖ʐς̍쐬
  mesh().createFaceAreas();

  // Face �̃��C�g�̍쐬
  if ( mesh().face_mates().empty() )
    mesh().createFaceMates();

  // Dual Graph �̍쐬
  dg_.setMesh( mesh() );
  dg_.createDG();

  // �k�ނ����ʂ̃m�[�h���폜����
  dg_.checkDG();

  // Face ���C�g�͂�������Ȃ�����폜
  mesh().clearFaceMates();

  cout << "initialize done. " << endl;
}

//
// �|���S���̖ʂ� Slimball �ɕϊ�
//
void SlimSimplify::polyToSlimBall()
{
  cout << "convert original polygons to slimballs ... " << endl;
  std::vector<DGNode*>& nodes = dg_.nodes();
  for ( int i = 0; i < nodes.size(); ++i )
    {
      DGNode* node = nodes[i];
      if ( node->isDeleted() ) continue;

      // SlimBall �̍쐬
      SlimBalld* slimball = createSlimBall();
      
      // DGNode �� SlimBall ��ǉ�
      node->setSlimBall( slimball );

      PIEMHandler piem_handle( mesh() );

      // �T�|�[�g�̒��S�Ɣ��a���v�Z
      Point3d center;
      double radius;
      std::vector<unsigned int>& face_list = node->face_list();
      piem_handle.calcCenterRadius( face_list, center, &radius );

      // �񎟊֐��̌W�����v�Z
      Quadraticd* quad = new Quadraticd;
      std::vector<double> coeff_global( QUADRATIC_COEFF );
      // Mesh �̖ʂ͈�����o�^����Ă��Ȃ��͂�
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

      // �o�^
      slimball->setCenter( center );
      slimball->setSupport( radius );
      slimball->addBasisFunction( *quad );

      // �덷��o�^���Ă���
      slimball->setUserDefined( 0.0 );
    }
  cout << "convert original polygons to slimballs done. " << endl;
}

double SlimSimplify::optimize( DGEdge* edge, bool isUpdate )
{
  if ( edge->slimball() ) edge->deleteSlimBall();
  if ( edge->piem() ) edge->deletePIEM();
  
  // SlimBall �̍쐬
  SlimBalld* slimball = new SlimBalld;
      
  // PIEM �̍쐬
  PIEM* piem = new PIEM;

  // ���[�̃m�[�h���� face list ���擾
  std::vector<unsigned int> face_list;
  DGNode* sn = edge->sn();
  std::vector<unsigned int>& sfl = sn->face_list();
  for ( int i = 0; i < sfl.size(); ++i ) { face_list.push_back( sfl[i] ); };
  DGNode* en = edge->en();
  std::vector<unsigned int>& efl = en->face_list();
  for ( int i = 0; i < efl.size(); ++i ) { face_list.push_back( efl[i] ); };

//       cout << "face_list " << face_list.size() << endl;

  // PIEM �ɂ��œK�֐��̌v�Z

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

  // crease �̂Ƃ��͌덷�֐��l�� 1000 �{�ɂ���
//   if ( edge->isCrease() == true ) error *= 1000.0;
  if ( edge->isCrease() == true ) error *= CREASE_WEIGHT;
      
  slimball->setCenter( center );
  slimball->setSupport( radius );
  slimball->addBasisFunction( *quad );
  // �덷��o�^���Ă���
  slimball->setUserDefined( error );

  // update �̒��Ŏg���鎞�ȊO�� slimball ����������
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

// ���_�̃m�[�h�� PIEM ���쐬
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

      // �m�[�h�� face_list �͂�������Ȃ��i�͂��j�Ȃ̂ŏ���
      node->deleteFaceList();
    }
  cout << "create node PIEMs done. " << endl;
}

void SlimSimplify::setPQ()
{
  cout << "create priority queue ... " << endl;

  // queue �̏�����
  pq_.init( dg_.edges().size() );

//   cout << " isStorePIEM " << isStorePIEM() << endl;
  if ( isStorePIEM() )
    {
      createNodePIEMs();
    }

  // MeshDG �̃G�b�W���� PQ �̃m�[�h���쐬
  for ( int k = 0; k < dg_.edges().size(); ++k )
    {
      DGEdge* edge = dg_.edge(k);
      if ( edge->isDeleted() ) continue;
      
      double error = optimize( edge, false );
      
//       cout << "edge " << k << " error " << error << endl;

      // PQ �m�[�h�̍쐬
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

  // slimball �����X�g�ɒǉ�
  addSlimBall( edge->slimball() );

//   edge->slimball()->print();

  // SlimBall �̐e�q�֌W�̍\�z
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

  // �G�b�W�� slimball �� PIEM �� sn �Ɉ����p��
  // sn �ɂ��Ƃ��ƕt���Ă��� PIEM �͏���, slimball �̓����N��؂�
  
  sn->setSlimBall( edge->slimball() );
  
  sn->deletePIEM();
  sn->setPIEM( edge->piem() );

  // ���p���I������̂� NULL �ɂ��Ă���
  edge->setSlimBall( NULL );
  edge->setPIEM( NULL );

  // en �� face_list �� edges �� sn �ɒǉ�
  sn->add( en, edge );

//   cout << "\t new node ... num of faces: " << sn->face_list().size() << endl;

  // en �� edge �� clear
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

      // �������ׂ��m�[�h�̃T�C�Y�� ���Ƃ̖ʂ̃T�C�Y�� 1/10 �ɂȂ�����
      // PIEM �ۑ��ɐ؂�ւ���
//       cout << "tmp " << slimballs_.size() - mesh().numFaces() << endl;
      if ( ( count > num_faces * .9 ) && ( changeNodeStoreType == false ) && ( isStorePIEM() == false) )
	{
	  createNodePIEMs();
	  setIsStorePIEM( true );
	  changeNodeStoreType = true;
	}

      DGNode* sn = edge->sn();

      update( edge );

      // �m�[�h�ɗאڂ���G�b�W�� slimball �� PIEM ���A�b�v�f�[�g
      std::vector<DGEdge*> edges = sn->edges();
      for ( int i = 0; i < edges.size(); ++i )
	{
	  DGEdge* ne = edges[i];
	  if ( ne->isDeleted() ) continue;

	  // slimball �� PIEM ���A�b�v�f�[�g
	  double error = optimize( ne, false );
	  PQNoded pqnode( ne->id(), error );
	  pq_.update( ne->id(), pqnode );
	}
    }

  cout << "simplification done. " << endl;
}

// slimball �̕��ёւ�
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
      // level 0 �� slimball ���擾
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
		  // �t�m�[�h��ۑ����Ȃ��ꍇ�͎��̍s���R�����g�A�E�g����
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
