#include "decl.h"
using namespace std;
using namespace Eigen;


int patchgen::df_check_geometry(GeometryData& gd)
{
	std::vector<std::vector<Eigen::VectorXi>> tmp_qg;

	for (int i = 0; i < gd.quadrangulable_geometry.size(); i++) {
		std::cout << "Edge_No. " << i << "\n";
		for (int j = 0; j < gd.quadrangulable_geometry[i].size(); j++) {
			for (int k = 0; k < gd.quadrangulable_geometry[i][j].size(); k++) {
				std::cout << gd.quadrangulable_geometry[i][j](k) << " ";
			}
			std::cout << "\n";
		}
	}

	for (int i = 0; i < gd.quadrangulable_geometry.size(); i++) {
		std::vector<Eigen::VectorXi> tmp_tmp_qg;
		for (int j = 0; j < gd.quadrangulable_geometry[i].size(); j++) {
			int check_s_o_e = df_check_singularity_on_edge(gd.subpatch_leaf[i], gd.quadrangulable_geometry[i][j]);
			if (check_s_o_e == 0)
				tmp_tmp_qg.push_back(gd.quadrangulable_geometry[i][j]);
		}
		tmp_qg.push_back(tmp_tmp_qg);
	}

	gd.quadrangulable_geometry = tmp_qg;

	std::cout << "\n";

	for (int i = 0; i < gd.quadrangulable_geometry.size(); i++) {
		std::cout << "Edge_No. " << i << "\n";
		for (int j = 0; j < gd.quadrangulable_geometry[i].size(); j++) {
			for (int k = 0; k < gd.quadrangulable_geometry[i][j].size(); k++) {
				std::cout << gd.quadrangulable_geometry[i][j](k) << " ";
			}
			std::cout << "\n";
		}
	}
	
	return df_check_quadrangulable_geometry(gd);
}

int patchgen::df_check_singularity_on_edge(std::vector<SubPatchData> s_leaf, Eigen::VectorXi variables)
{
	/*
	std::cout << "variables = ";
	for (int k = 0; k < variables.size(); k++) {
		std::cout << variables(k) << " ";
	}
	std::cout << "\n";
	*/
	for (int i = 0; i < s_leaf.size(); i++) {

		int num_sub_sides = s_leaf[i].num_sides;

		Eigen::VectorXi X;
		X.resize(num_sub_sides);

		for (int j = 0; j < s_leaf[i].side_subdivi.size(); j++) {
			if (s_leaf[i].side_subdivi[j]->num_subdivi == -1) {
				if (s_leaf[i].side_subdivi[j]->ilp_variable_idx != -1) {
					X[j] = variables(s_leaf[i].side_subdivi[j]->ilp_variable_idx);
				}
				else {
					SubPatchEdgeData *t;
					t = s_leaf[i].side_subdivi[j]->same_edge_pointer;
					for (;;) {
						if (t->same_edge_pointer == NULL) break;
						t = t->same_edge_pointer;
					}
					if (t->num_subdivi != -1)
						X[j] = t->num_subdivi;
					else
						X[j] = variables(t->ilp_variable_idx);
				}
			}
			else {
				X[j] = s_leaf[i].side_subdivi[j]->num_subdivi;
			}
		}
		//std::cout << "s_leaf idx = " << s_leaf[i].idx << "\n";
		if (num_sub_sides == 5) {
			Eigen::MatrixXi Ai;
			Ai.resize(num_sub_sides, num_sub_sides);

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
			/*
			std::cout << "s_leaf(subdivi) = ";
			for (int k = 0; k < X.size(); k++) {
				std::cout << X(k) << " ";
			}
			std::cout << "\n";
			*/
			X = Ai * X;
			X = X / 2;
		}
		/*
		std::cout << "s_leaf(polychord) = ";
		for (int k = 0; k < X.size(); k++) {
			std::cout << X(k) << " ";
		}
		std::cout << "\n";
		*/
		for (int j = 0; j < X.size(); j++) 
			if (X[j] == 0) return 1;
	}

	return 0;
}

int patchgen::df_check_quadrangulable_geometry(GeometryData& gd)
{
	for (int i = 0; i < gd.quadrangulable_geometry.size(); i++) 
		if (gd.quadrangulable_geometry[i].size() != 0) 
			return 0;
	return 1;
}

int patchgen::df_check_ilp_side_subdivi(Eigen::VectorXi ilp_side_subdivi)
{
	for (int i = 0; i < ilp_side_subdivi.size(); i++)
		if (ilp_side_subdivi[i] == 0)
			return 1;
	return 0;
}

