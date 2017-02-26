#include "decl.h"
#include <kt84/util.h>
#include <kt84/openmesh/flip_faces.hh>
#include "../patchgen_demo/curvenetwork/Patch.hh"
#include "../patchgen_demo/curvenetwork/decl.hh"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <kt84/eigen_def.hh>
#include <kt84/openmesh/edgeloop.hh>
#include <sstream>


int patchgen::df_simple_pentagon(GeometryData& gd)
{
	std::cout << "-----Simple Pentagon-----\n";

	int num_variables = gd.num_sides; //ポリコードの数

	Eigen::MatrixXi Ai;
	Eigen::VectorXi X;
	Ai.resize(num_variables, num_variables);
	X.resize(num_variables, num_variables);
	
	//A << 1, 1, 0, 0, 0,
	//	   0, 0, 1, 1, 0,
	//	   1, 0, 0, 0, 1,
	//	   0, 1, 0, 1, 0,
	//	   0, 0, 1, 0, 1;

	Ai << 1, 1, 1, -1, -1,
		-1, 1, 1, 1, -1,
		-1, -1, 1, 1, 1,
		1, -1, -1, 1, 1,
		1, 1, -1, -1, 1;

	gd.side_subdivi = gd.ilp_side_subdivi;

	X = Ai * gd.side_subdivi;
	X = X / 2;

	gd.patch_geometry = X;


	return 0;
}

int patchgen::df_set_simplepatch(GeometryData& gd)
{
	gd.subpatch.clear();

	gd.side_subdivi = gd.ilp_side_subdivi;
	SubPatchData *s;
	s = new SubPatchData;
	IdxCounter ic;
	std::vector<SubPatchEdgeData*> s_s;
	s->num_sides = gd.num_sides;
	s->c_pointer_r = s->c_pointer_l = s->p_pointer = NULL;
	s->idx = ic.getPatchIdx();
	
	for (int i = 0; i < gd.num_sides; i++) {
		SubPatchEdgeData *se;
		se = new SubPatchEdgeData;
		se->belong_patch = s;
		se->c_pointer_r = se->c_pointer_l = se->p_pointer = se->same_edge_pointer = NULL;
		se->idx = ic.getEdgeIdx();
		se->num_subdivi = gd.side_subdivi[i];
		s_s.push_back(se);
	}	
	s->side_subdivi = s_s;
	df_debug_diplay_SubPatchData(*s);///////
	gd.subpatch.push_back(*s);
	
	//5角形でも特異点が複数必要な場合は調整必要
	std::vector<SubPatchData> tmp_leaf;
	tmp_leaf.push_back(*s);
	gd.subpatch_leaf.push_back(tmp_leaf);
	//

	gd.quadrangulable_geometry.clear();
	gd.quadrangulable_geometry.resize(1);
	gd.quadrangulable_geometry[0].resize(1);
	gd.quadrangulable_geometry[0][0] = gd.side_subdivi;

	return 0;
}