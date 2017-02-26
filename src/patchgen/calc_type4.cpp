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
#include "ILP.h"


int patchgen::df_over_hexagon(GeometryData& gd)
{
	std::cout << "-----Over Hexagon-----\n";

	int num_sides = gd.num_sides;
	int num_oppside = (int)floor(num_sides / 2);

	gd.subpatch.clear();
	gd.subpatch_leaf.clear();
	gd.subpatch_leaf.resize(num_sides);
	gd.quadrangulable_geometry.clear();
	gd.quadrangulable_geometry.resize(num_sides);

	//最初に分割する対辺を指定
	for (int i = 0; i < num_sides; i++) {
		IdxCounter ic;
		std::vector<SubPatchData*> s_leaf;
		SubPatchData *s_root;
		s_root = new SubPatchData;
		s_root->num_sides = gd.num_sides;
		for (int i = 0; i < gd.side_subdivi.size(); i++) {
			SubPatchEdgeData *se;
			se = new SubPatchEdgeData;
			se->num_subdivi = gd.side_subdivi[i];
			se->c_pointer_l = se->c_pointer_r = se->p_pointer = se->same_edge_pointer = NULL;
			se->belong_patch = s_root;
			se->idx = ic.getEdgeIdx();
			s_root->side_subdivi.push_back(se);
		}
		s_root->c_pointer_l = s_root->c_pointer_r = s_root->p_pointer = NULL;
		s_root->idx = ic.getPatchIdx();
		s_leaf.push_back(s_root);

		//std::cout << "-----Set SubPatchData-----\n";
		df_set_subpatch(s_leaf, ic, i);

		OverHexagonILP oc;
		oc.num_variables = (gd.num_sides - 5) * 3;
		df_set_topology_over_hexagon(*s_root, gd, oc);
		df_solve_topology_over_hexagon(gd, oc, i);

		//df_debug_diplay_SubPatchData(*s_root);

		//deleteしてね！！！！！！！！!！！！！！！！！！！！！！！！
		//メモリリークしちゃうよ！！！！！！！！！！！！！！！！！！！
		gd.subpatch_leaf[i] = oc.s_leaf;
		gd.subpatch.push_back(*s_root);
	}
	
	/*
	for (int i = 0; i < gd.subpatch_leaf.size(); i++) {
		std::cout << "Leaf_Edge_No. " << i << "\n";
		for (int j = 0; j < gd.subpatch_leaf[i].size(); j++) {
			std::cout << "leaf_idx = " << gd.subpatch_leaf[i][j].idx << "\n";
		}
	}
	*/
	/*
	for (int i = 0; i < gd.quadrangulable_geometry.size(); i++) {
		std::cout << "Edge_No. " << i <<"\n";
		for (int j = 0; j < gd.quadrangulable_geometry[i].size(); j++) {
			for (int k = 0; k < gd.quadrangulable_geometry[i][j].size(); k++) {
				std::cout << gd.quadrangulable_geometry[i][j](k) << " ";
			}
			std::cout << "\n";
		}
	}
	*/
	std::cout << "----------fin----------\n";
	return 0;
}

Eigen::VectorXi patchgen::df_calc_ilp_pentagon_polychord(SubPatchData* s, GeometryData& gd)
{
	//quadrangulable_geometryが空のときはエラー起きる
	int num_sub_sides = s->side_subdivi.size();
		
	Eigen::VectorXi ip;
	ip.resize(num_sub_sides);

	for (int i = 0; i < num_sub_sides; i++) {
		if (s->side_subdivi[i]->num_subdivi == -1) {
			if (s->side_subdivi[i]->ilp_variable_idx != -1) {
				s->side_subdivi[i]->ilp_num_subdivi = gd.quadrangulable_geometry[gd.subdivi_edge_no][gd.geometry_no](s->side_subdivi[i]->ilp_variable_idx);
			}
			else {
				SubPatchEdgeData *t;
				t = s->side_subdivi[i]->same_edge_pointer;
				for (;;) {
					if (t->same_edge_pointer == NULL) break;
					t = t->same_edge_pointer;
				}
				if (t->num_subdivi != -1)
					s->side_subdivi[i]->ilp_num_subdivi = t->num_subdivi;
				else
					s->side_subdivi[i]->ilp_num_subdivi = gd.quadrangulable_geometry[gd.subdivi_edge_no][gd.geometry_no](t->ilp_variable_idx);
			}
		}
		else {
			s->side_subdivi[i]->ilp_num_subdivi = s->side_subdivi[i]->num_subdivi;
		}
		ip(i) = s->side_subdivi[i]->ilp_num_subdivi;
	}

	if (num_sub_sides == 5) {
		Eigen::MatrixXi Ai;
		Ai.resize(num_sub_sides, num_sub_sides);

		//A << 1, 1, 0, 0, 0,
		//	   0, 0, 1, 1, 0,
		//	   1, 0, 0, 0, 1,
		//	   0, 1, 0, 1, 0,
		//	   0, 0, 1, 0, 1;

		Ai << 1,  1,  1, -1, -1,
			 -1,  1,  1,  1, -1,
			 -1, -1,  1,  1,  1,
			  1, -1, -1,  1,  1,
		 	  1,  1, -1, -1, 1;
	
		ip = Ai * ip;
		ip = ip / 2;
	}

	return ip;
}


int patchgen::df_set_subpatch(std::vector<SubPatchData*>& s_node, IdxCounter& ic, int start_edge)
{
	std::cout << "-----Set SubPatchData-----\n";

	int num_width = s_node.size();
	bool all_pentagon_flag = true;
	std::vector<SubPatchData*> s_leaf;
	for (int node_count = 0; node_count < num_width; node_count++) {

		int num_sides = s_node[node_count]->num_sides;

		if (num_sides > 5) {
			int num_oppside = (int)floor(num_sides / 2);

			SubPatchData *s_r;
			s_r = new SubPatchData;
			s_r->num_sides = num_oppside + 2;
			s_r->idx = ic.getPatchIdx();
			SubPatchData *s_l;
			s_l = new SubPatchData;
			s_l->num_sides = num_sides - num_oppside + 2;
			s_l->idx = ic.getPatchIdx();

			//分割するエッジの開始インデックス
			int sub_line = num_sides - 1;
			int side_count = (s_node[node_count]->p_pointer == NULL) ? start_edge : num_sides - 1;////

			for (int sc = 0; sc < s_r->num_sides - 1; sc++, side_count++) {
				SubPatchEdgeData *se;
				se = new SubPatchEdgeData;
				if (sc == 0) {
					se->num_subdivi = -1;
					se->p_pointer = s_node[node_count]->side_subdivi[side_count];
					s_node[node_count]->side_subdivi[side_count]->c_pointer_r = se;
					se->same_edge_pointer = NULL;
				} else if(sc == s_r->num_sides - 2) {
					se->num_subdivi = -1;
					se->p_pointer = s_node[node_count]->side_subdivi[side_count];
					s_node[node_count]->side_subdivi[side_count]->c_pointer_r = se;
					se->same_edge_pointer = NULL;
				} else {
					se->num_subdivi = s_node[node_count]->side_subdivi[side_count]->num_subdivi;
					se->p_pointer = NULL;
					se->same_edge_pointer = s_node[node_count]->side_subdivi[side_count];
				}
				se->belong_patch = s_r;
				se->idx = ic.getEdgeIdx();
				se->c_pointer_r = se->c_pointer_l = NULL;
				s_r->side_subdivi.push_back(se);
				side_count = (side_count == num_sides - 1) ? -1 : side_count;
			}

			SubPatchEdgeData *se_r;
			se_r = new SubPatchEdgeData;
			se_r->num_subdivi = -1;
			se_r->p_pointer = se_r->c_pointer_r = se_r->c_pointer_l = se_r->same_edge_pointer = NULL;
			se_r->belong_patch = s_r;
			se_r->idx = ic.getEdgeIdx();

			s_r->side_subdivi.push_back(se_r);
			side_count = (side_count == 0) ? num_sides : side_count;
			side_count--;

			for (int sc = 0; sc < s_l->num_sides - 1; sc++, side_count++) {
				SubPatchEdgeData *se;
				se = new SubPatchEdgeData;
				if (sc == 0) {
					se->num_subdivi = -1;
					se->p_pointer = s_node[node_count]->side_subdivi[side_count];
					s_node[node_count]->side_subdivi[side_count]->c_pointer_l = se;
					se->same_edge_pointer = NULL;
				}
				else if (sc == s_l->num_sides - 2) {
					se->num_subdivi = -1;
					se->p_pointer = s_node[node_count]->side_subdivi[side_count];
					s_node[node_count]->side_subdivi[side_count]->c_pointer_l = se;
					se->same_edge_pointer = NULL;
				}
				else {
					se->num_subdivi = s_node[node_count]->side_subdivi[side_count]->num_subdivi;
					se->p_pointer = NULL;
					se->same_edge_pointer = s_node[node_count]->side_subdivi[side_count];
				}
				se->belong_patch = s_l;
				se->idx = ic.getEdgeIdx();
				se->c_pointer_r = se->c_pointer_l = NULL;
				s_l->side_subdivi.push_back(se);
				side_count = (side_count == num_sides - 1) ? -1 : side_count;
			}

			SubPatchEdgeData *se_l;
			se_l = new SubPatchEdgeData;
			se_l->num_subdivi = -1;
			se_l->p_pointer = se_l->c_pointer_r = se_l->c_pointer_l = se_l->same_edge_pointer = NULL;
			se_l->belong_patch = s_l;
			se_l->idx = ic.getEdgeIdx();

			s_l->side_subdivi.push_back(se_l);

			s_node[node_count]->c_pointer_r = s_r;
			s_node[node_count]->c_pointer_l = s_l;
			s_r->p_pointer = s_l->p_pointer = s_node[node_count];
			s_r->c_pointer_r = s_r->c_pointer_l = s_l->c_pointer_r = s_l->c_pointer_l = NULL;
			
			s_leaf.push_back(s_r);
			s_leaf.push_back(s_l);

			all_pentagon_flag = false;
		}
	}

	return all_pentagon_flag ? 0 : df_set_subpatch(s_leaf, ic, start_edge);
}

void patchgen::df_debug_diplay_SubPatchData(SubPatchData s)
{
	std::cout << "SubPatch_idx : " << s.idx << "\n";
	std::cout << "		patch_parent_idx : " << ((s.p_pointer != NULL) ? s.p_pointer->idx : -1) << "\n";
	std::cout << "		patch_childR_idx : " << ((s.c_pointer_r != NULL) ? s.c_pointer_r->idx : -1) << "\n";
	std::cout << "		patch_childL_idx : " << ((s.c_pointer_l != NULL) ? s.c_pointer_l->idx : -1) << "\n";
	
	std::cout << "	num_subdivi : " << s.num_sides << "\n";
	
	for (int i = 0; i < s.side_subdivi.size(); i++) {
		
		std::cout << "		Edge_idx : " << s.side_subdivi[i]->idx << "		" << s.side_subdivi[i]->num_subdivi << "\n";
		std::cout << "			edge_parent_idx  : " << ((s.side_subdivi[i]->p_pointer != NULL) ? s.side_subdivi[i]->p_pointer->idx : -1) << "\n";
		std::cout << "			edge_childR_idx  : " << ((s.side_subdivi[i]->c_pointer_r != NULL) ? s.side_subdivi[i]->c_pointer_r->idx : -1) << "\n";
		std::cout << "			edge_childL_idx  : " << ((s.side_subdivi[i]->c_pointer_l != NULL) ? s.side_subdivi[i]->c_pointer_l->idx : -1) << "\n";
		std::cout << "			edge_same_idx    : " << ((s.side_subdivi[i]->same_edge_pointer != NULL) ? s.side_subdivi[i]->same_edge_pointer->idx : -1) << "\n";
		std::cout << "			belong_patch_idx : " << ((s.side_subdivi[i]->belong_patch != NULL) ? s.side_subdivi[i]->belong_patch->idx : -1) << "\n";
		std::cout << "			ilp_variable_idx : " << s.side_subdivi[i]->ilp_variable_idx << "\n";
	}
	
	
	if (s.c_pointer_r == NULL)
		;
	else
		df_debug_diplay_SubPatchData(*s.c_pointer_r);

	if (s.c_pointer_l == NULL)
		;
	else
		df_debug_diplay_SubPatchData(*s.c_pointer_l);

}


int patchgen::df_set_topology_over_hexagon(SubPatchData& s_root, GeometryData& gd, OverHexagonILP& oc)
{
	std::cout << "-----Set topology for over hexagon-----\n";

	int num_subpatch = 3 + 2 * (gd.num_sides - 6) - 1;//subpatchの総数は 3 + 2 * (num_sides - 6)からroot分を引く，rootは自明なので計算しない
	int num_edge_variables = num_subpatch * 3;
	int num_variables = num_subpatch + num_edge_variables;
	std::vector<SubPatchData> s_leaf;
	std::vector<Eigen::Vector2i> s_consis;

	//----------------loopして範囲内の分割数を列挙------------------------
	//for

	//4_各サブパッチはエッジの合計が偶数である 20160908->まだ
	std::vector<SubPatchData> s_now;
	s_now.push_back(*s_root.c_pointer_r);
	s_now.push_back(*s_root.c_pointer_l);

	Eigen::MatrixXd cons4;
	Eigen::VectorXd cons4_rhs;
	cons4.resize(num_subpatch, num_variables);
	cons4_rhs.resize(num_subpatch);

	for (int i = 0; i < num_subpatch; i++)
		for (int j = 0; j < num_variables; j++)
			cons4(i, j) = 0;
	for (int i = 0; i < num_subpatch; i++)
		cons4_rhs(i) = 0;
	//サブパッチをひとつ増やすごとに新たに追加される変数は3つ
	int subpatch_count = 0;
	for (;;) {
		std::vector<SubPatchData> s_now_tmp;
		for (int i = 0; i < s_now.size(); i++) {
			int sp_variables = 0;
			for (int j = 0; j < s_now[i].side_subdivi.size(); j++) {
				if (s_now[i].side_subdivi[j]->num_subdivi == -1 && s_now[i].side_subdivi[j]->same_edge_pointer == NULL) {
					s_now[i].side_subdivi[j]->ilp_variable_idx = (s_now[i].idx - 1) * 3 + sp_variables;
					cons4(s_now[i].idx - 1, s_now[i].side_subdivi[j]->ilp_variable_idx) = 1;
					sp_variables++;
				}
				else if (s_now[i].side_subdivi[j]->num_subdivi == -1 && s_now[i].side_subdivi[j]->same_edge_pointer != NULL) {
					SubPatchEdgeData *t;
					t = s_now[i].side_subdivi[j]->same_edge_pointer;
					for (;;) {
						if (t->same_edge_pointer == NULL) break;
						t = t->same_edge_pointer;
					}
					if (t->num_subdivi == -1)
						cons4(s_now[i].idx - 1, t->ilp_variable_idx) = 1;
					else
						cons4_rhs(s_now[i].idx - 1) -= t->num_subdivi;
				}
				else {
					cons4_rhs(s_now[i].idx - 1) -= s_now[i].side_subdivi[j]->num_subdivi;
				}
			}
			cons4(s_now[i].idx - 1, num_subpatch * 3 + subpatch_count++) = -2;
			if (s_now[i].c_pointer_r != NULL) s_now_tmp.push_back(*s_now[i].c_pointer_r);
			if (s_now[i].c_pointer_l != NULL) s_now_tmp.push_back(*s_now[i].c_pointer_l);
				
			int new_flag = -1;
			if (s_consis.size() == 0)
				new_flag = -1;
			else
				for (int s = 0; s < s_consis.size(); s++)
					if (s_consis[s](1) == s_now[i].side_subdivi.size())
						new_flag = s;

			if (new_flag == -1) {
				Eigen::Vector2i t;
				t << s_now[i].side_subdivi.size(), 0;
				s_consis.push_back(t);
			} else {
				s_consis[new_flag](1)++;
			}

			if (s_now[i].c_pointer_r == NULL && s_now[i].c_pointer_l == NULL) s_leaf.push_back(s_now[i]);
		}
		if (s_now_tmp.size() == 0) break;
		else s_now = s_now_tmp;
	}
	/*
	std::cout << "---4------sahen------\n";
	std::cout << cons4 << "\n";
	std::cout << "---4------uhen------\n";
	std::cout << cons4_rhs << "\n";
	*/


	//0_シンプルである

	int num_simple = 0;
	for (auto c : s_consis) {
		num_simple += c(0);
		Eigen::Vector2i t;
		t << 0, 12;
		oc.ilp_counter_slack.push_back(t);
	}

	s_now.clear();
	s_now.push_back(*s_root.c_pointer_r);
	s_now.push_back(*s_root.c_pointer_l);

	Eigen::MatrixXd cons0;
	Eigen::VectorXd cons0_rhs;
	cons0.resize(num_simple, num_variables);
	cons0_rhs.resize(num_simple);
	Eigen::VectorXd cons0_ns_rhs;
	cons0_ns_rhs.resize(num_simple);

	for (int i = 0; i < num_simple; i++)
		for (int j = 0; j < num_variables; j++)
			cons0(i, j) = 0;

	for (int i = 0; i < num_simple; i++)
		cons0_rhs(i) = cons0_ns_rhs(i) = 0;

	int tate = 0;
	for (;;) {
		std::vector<SubPatchData> s_now_tmp;
		for (int i = 0; i < s_now.size(); i++) {
			int num_subsides = s_now[i].side_subdivi.size();
			for (int s = 0; s < num_subsides; s++) {
				cons0_rhs(s + tate) = (num_subsides == 5) ? 1 : (num_subsides - 4) * 2;
				cons0_ns_rhs(s + tate) = (num_subsides == 5) ? 0 : (num_subsides - 4) * 2;
			}
			for (int j = 0; j < num_subsides; j++) {
				for (int k = 0; k < num_subsides; k++) {
					int c = (j + k >= num_subsides) ? (j + k - num_subsides) : (j + k);
					if (s_now[i].side_subdivi[c]->num_subdivi == -1 && s_now[i].side_subdivi[c]->same_edge_pointer == NULL) {
						cons0(tate, s_now[i].side_subdivi[c]->ilp_variable_idx) = (k < (num_subsides - 2)) ? 1 : -1;
					}
					else if (s_now[i].side_subdivi[c]->num_subdivi == -1 && s_now[i].side_subdivi[c]->same_edge_pointer != NULL) {
						SubPatchEdgeData *t;
						t = s_now[i].side_subdivi[c]->same_edge_pointer;
						for (;;) {
							if (t->same_edge_pointer == NULL) break;
							t = t->same_edge_pointer;
						}
						if (t->num_subdivi == -1) {
							cons0(tate, t->ilp_variable_idx) = (k < (num_subsides - 2)) ? 1 : -1;
						} else {
							cons0_rhs(tate) += (k < (num_subsides - 2)) ? -t->num_subdivi : t->num_subdivi;
							cons0_ns_rhs(tate) += (k < (num_subsides - 2)) ? -t->num_subdivi : t->num_subdivi;
						}
					}
					else {
						cons0_rhs(tate) += (k < (num_subsides - 2)) ? -s_now[i].side_subdivi[c]->num_subdivi : s_now[i].side_subdivi[c]->num_subdivi;
						cons0_ns_rhs(tate) += (k < (num_subsides - 2)) ? -s_now[i].side_subdivi[c]->num_subdivi : s_now[i].side_subdivi[c]->num_subdivi;
					}
				}
				tate++;
			}
			if (s_now[i].c_pointer_r != NULL) s_now_tmp.push_back(*s_now[i].c_pointer_r);
			if (s_now[i].c_pointer_l != NULL) s_now_tmp.push_back(*s_now[i].c_pointer_l);
		}
		if (s_now_tmp.size() == 0) break;
		else s_now = s_now_tmp;
	}

	//5_分割する場所の条件
	int num_bunkatsu = (gd.num_sides - 5) * 3;
	std::vector<SubPatchData> s_now2;
	s_now2.push_back(s_root);

	Eigen::MatrixXd cons5;
	Eigen::VectorXd cons5_rhs;
	cons5.resize(num_bunkatsu, num_variables);
	cons5_rhs.resize(num_bunkatsu);

	for (int i = 0; i < num_bunkatsu; i++)
		for (int j = 0; j < num_variables; j++)
			cons5(i, j) = 0;
	for (int i = 0; i < num_bunkatsu; i++)
		cons5_rhs(i) = 0;

	int depth = 0;//いらない
	int depth2 = 0;
	for (;;) {
		std::vector<SubPatchData> s_now_tmp;
		for (int i = 0; i < s_now2.size(); i++) {
			int bunkatsu_count = 0;
			for (int j = 0; j < s_now2[i].side_subdivi.size(); j++){
				if (s_now2[i].side_subdivi[j]->c_pointer_r != NULL && s_now2[i].side_subdivi[j]->c_pointer_l != NULL) {

					if (s_now2[i].side_subdivi[j]->num_subdivi == -1)
						cons5(3 * depth2 + bunkatsu_count, s_now2[i].side_subdivi[j]->ilp_variable_idx) = -1;
					else
						cons5_rhs(3 * depth2 + bunkatsu_count) += s_now2[i].side_subdivi[j]->num_subdivi;

					cons5(3 * depth2 + bunkatsu_count, s_now2[i].side_subdivi[j]->c_pointer_r->ilp_variable_idx) = 1;
					cons5(3 * depth2 + bunkatsu_count, s_now2[i].side_subdivi[j]->c_pointer_l->ilp_variable_idx) = 1;
					bunkatsu_count++;

					//----------------
					Eigen::Vector2i t;
					t(0) = 1;
					t(1) = (s_now2[i].side_subdivi[j]->num_subdivi == -1) ? 10 : s_now2[i].side_subdivi[j]->num_subdivi;////////////////
					oc.ilp_counter.push_back(t);
					//----------------


					//----------------
					if (depth == 0) oc.selected_subdivi_edge.push_back(s_now2[i].side_subdivi[j]->num_subdivi);
					//----------------
				}
				if (bunkatsu_count == 2) break;
			}
			
			for (int k = 0; k < s_now2[i].c_pointer_r->side_subdivi.size(); k++) {
				if (s_now2[i].c_pointer_r->side_subdivi[k]->p_pointer == NULL && s_now2[i].c_pointer_r->side_subdivi[k]->num_subdivi == -1 && s_now2[i].c_pointer_r->side_subdivi[k]->same_edge_pointer == NULL) {
					cons5(3 * depth2 + bunkatsu_count, s_now2[i].c_pointer_r->side_subdivi[k]->ilp_variable_idx) = 1;

					//--------------------
					if (depth == 0) oc.first_subdivi_edge = s_now2[i].c_pointer_r->side_subdivi[k]->ilp_variable_idx;
					//--------------------

					//--------------------
					Eigen::Vector2i t;
					t(0) = 1;
					t(1) = 10;//////////////////////
					oc.ilp_counter.push_back(t);
					//--------------------
				}
			}

			for (int k = 0; k < s_now2[i].c_pointer_l->side_subdivi.size(); k++) {
				if (s_now2[i].c_pointer_l->side_subdivi[k]->p_pointer == NULL && s_now2[i].c_pointer_l->side_subdivi[k]->num_subdivi == -1 && s_now2[i].c_pointer_l->side_subdivi[k]->same_edge_pointer == NULL) 
					cons5(3 * depth2 + bunkatsu_count, s_now2[i].c_pointer_l->side_subdivi[k]->ilp_variable_idx) = -1;
			}
			
			if (s_now2[i].c_pointer_r->c_pointer_r != NULL && s_now2[i].c_pointer_r->c_pointer_l != NULL) s_now_tmp.push_back(*s_now2[i].c_pointer_r);
			if (s_now2[i].c_pointer_l->c_pointer_r != NULL && s_now2[i].c_pointer_l->c_pointer_l != NULL) s_now_tmp.push_back(*s_now2[i].c_pointer_l);

			depth2++;
		}
		depth++;
		if (s_now_tmp.size() == 0) break;
		else s_now2 = s_now_tmp;
	}
	/*
	std::cout << "---5------sahen------\n";
	std::cout << cons5 << "\n";
	std::cout << "---5------uhen------\n";
	std::cout << cons5_rhs << "\n";
	*/



	Eigen::MatrixXd cons3;
	Eigen::VectorXd cons3_rhs;
	cons3.resize(num_edge_variables, num_variables);
	cons3_rhs.resize(num_edge_variables);

	for (int i = 0; i < num_edge_variables; i++)
		for (int j = 0; j < num_variables; j++)
			cons3(i, j) = 0;
	for (int i = 0; i < num_edge_variables; i++)
		cons3_rhs(i) = 1;

	for (int i = 0; i < num_edge_variables; i++) {
		cons3(i, i) = 1;
		if ((i % 3) == 2) cons3_rhs(i) = 1;/*分割によって発生したエッジの最小分割数を2にしました*/
	}
	/*
	std::cout << "---3------sahen------\n";
	std::cout << cons3 << "\n";
	std::cout << "---3------uhen------\n";
	std::cout << cons3_rhs << "\n";
	*/
	oc.s_leaf = s_leaf;
	oc.cons0 = cons0;
	oc.cons0_rhs = cons0_rhs;
	oc.cons5 = cons5;
	oc.cons5_rhs = cons5_rhs;
	oc.cons4 = cons4;
	oc.cons4_rhs = cons4_rhs;
	oc.cons3 = cons3;
	oc.cons3_rhs = cons3_rhs;
	oc.cons0_ns_rhs = cons0_ns_rhs;

	return 0;
}



int patchgen::df_solve_topology_over_hexagon(GeometryData& gd, OverHexagonILP& oc, int selected_subdivi_edge)
{
	int num_subpatch = 3 + 2 * (gd.num_sides - 6) - 1;//subpatchの総数は 3 + 2 * (num_sides - 6)からroot分を引く，rootは自明なので計算しない
	int num_edge_variables = num_subpatch * 3;
	int num_variables = num_subpatch + num_edge_variables;

	Eigen::VectorXd cons_ex0;
	cons_ex0.resize(num_variables);
	for (int i = 0; i < num_variables; i++)
		cons_ex0(i) = (i == 0) ? 1 : 0;
	Eigen::VectorXd cons_ex1;
	cons_ex1.resize(num_variables);
	for (int i = 0; i < num_variables; i++)
		cons_ex1(i) = (i == 1) ? 1 : 0;
	Eigen::VectorXd cons_ex2;
	cons_ex2.resize(num_variables);
	for (int i = 0; i < num_variables; i++)
		cons_ex2(i) = (i == oc.first_subdivi_edge) ? 1 : 0;

	//オブジェクト関数のセット
	Eigen::VectorXd object;
	object.resize(num_variables);
	for (int i = 0; i < num_variables; i++) {
		if (i == oc.first_subdivi_edge)
			object(i) = 1;
		else
			object(i) = 0;
	}

	int max_subdivi_edge = 100;
	for (;;) {
		//制約代入
		ILP ilp(num_variables);
		ilp.add_constraint(oc.cons4, EQ, oc.cons4_rhs);
		ilp.add_constraint(oc.cons0, GE, oc.cons0_rhs);
		ilp.add_constraint(oc.cons5, EQ, oc.cons5_rhs);
		ilp.add_constraint(oc.cons3, GE, oc.cons3_rhs);
		ilp.add_constraint(cons_ex2, LE, max_subdivi_edge);
		ilp.set_objective(object, true);

		bool solve_flag = true;
		if (!ilp.solve()) solve_flag = false; //error処理してね！！！！！！！！！！！！

		if (!solve_flag) break;

		auto ilp_variables = ilp.get_variables();
		max_subdivi_edge = ilp_variables(oc.first_subdivi_edge) - 1;
		/*
		for (int i = 0; i < ilp_variables.size(); ++i)
			std::cout << ilp_variables[i] << " ";
		std::cout << "\n";
		*/
		df_check_repetition_geometry_data(gd.quadrangulable_geometry[selected_subdivi_edge], ilp_variables);

	}
	/*
	for (int subdivi_edge_end = 1; subdivi_edge_end < oc.selected_subdivi_edge[1]; subdivi_edge_end++) {
		for (int subdivi_edge_start = 1; subdivi_edge_start < oc.selected_subdivi_edge[0]; subdivi_edge_start++) {
			max_subdivi_edge = 100;
			for (;;) {
				//制約代入
				ILP ilp(num_variables);
				ilp.add_constraint(oc.cons4, EQ, oc.cons4_rhs);
				ilp.add_constraint(oc.cons0, GE, oc.cons0_ns_rhs);
				ilp.add_constraint(oc.cons5, EQ, oc.cons5_rhs);
				ilp.add_constraint(oc.cons3, GE, oc.cons3_rhs);
				ilp.add_constraint(cons_ex0, GE, subdivi_edge_start);
				ilp.add_constraint(cons_ex1, GE, subdivi_edge_end);
				ilp.add_constraint(cons_ex2, LE, max_subdivi_edge);
				ilp.set_objective(object, true);

				bool solve_flag = true;
				if (!ilp.solve()) solve_flag = false; //error処理してね！！！！！！！！！！！！
				auto ilp_variables = ilp.get_variables();

				if (!solve_flag) break;

				max_subdivi_edge = ilp_variables(oc.first_subdivi_edge) - 1;
				
				//for (int i = 0; i < ilp_variables.size(); ++i)
					//std::cout << ilp_variables[i] << " ";
				//std::cout << "\n";
				
				df_check_repetition_geometry_data(gd.quadrangulable_geometry[selected_subdivi_edge], ilp_variables);
			}
		}
	}
	*/
	return 0;
}


int patchgen::df_check_repetition_geometry_data(std::vector<Eigen::VectorXi>& quadrangulable_geometry, Eigen::VectorXi ilp_variables)
{
	if (quadrangulable_geometry.size() == 0) {
		quadrangulable_geometry.push_back(ilp_variables);
	} else {
		int check_num_i = 0;
		int num_i = quadrangulable_geometry.size();
		for (int i = 0; i < num_i; i++) {
			int check_num_j = 0;
			int num_j = quadrangulable_geometry[i].size();
			for (int j = 0; j < num_j; j++) {
				if (quadrangulable_geometry[i][j] == ilp_variables[j]) check_num_j++;
				else break;
			}
			if (num_j == check_num_j) break;
			else check_num_i++;
		}
		if (check_num_i == num_i) quadrangulable_geometry.push_back(ilp_variables);
	}
	return 0;
}