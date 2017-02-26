#pragma once
#include "PatchParam.h"
#include "GeometryData.h"
#include "SubPatchData.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <utility>


namespace patchgen {
	
	//ILPを計算 -> calc_ilp_sides.cpp
	Eigen::VectorXi df_calc_ilp_sides(Eigen::VectorXi& l);
	void df_debug_calc_ilp(std::vector<Eigen::VectorXi> l);

	//Polychord, パッチ分解の構造を計算 -> calc_quadrangulation.cpp
	void df_calc_quadrangulation(GeometryData& gd);
	int df_check_quadrangulation_type(GeometryData& gd);

	//シンプルなパッチを計算(345) -> calc_type1.cpp
	int df_simple_pentagon(GeometryData& gd);
	int df_set_simplepatch(GeometryData& gd);

	//シンプルなパッチを計算(345) -> calc_type1.cpp
	int df_over_hexagon(GeometryData& gd);
	int df_set_subpatch(std::vector<SubPatchData*>& s_node, IdxCounter& ic, int start_edge);
	void df_debug_diplay_SubPatchData(SubPatchData s);//debug!!!!!!!!!!!!!!!!!
	int df_set_topology_over_hexagon(SubPatchData& s, GeometryData& gd, OverHexagonILP& oc);
	int df_solve_topology_over_hexagon(GeometryData& gd, OverHexagonILP& oc, int selected_subdivi_edge);
	int df_check_repetition_geometry_data(std::vector<Eigen::VectorXi>& quadrangulable_geometry, Eigen::VectorXi ilp_variables);
	Eigen::VectorXi df_calc_ilp_pentagon_polychord(SubPatchData* s_node, GeometryData& gd);

	//sideの位置固定
	Eigen::Vector2d df_get_boundary_geometry(int num_sides, long double t, std::vector<Eigen::Vector2d> df_patch);

	//
	int df_check_geometry(GeometryData& gd);
	int df_check_singularity_on_edge(std::vector<SubPatchData> s_leaf, Eigen::VectorXi variables);
	int df_check_quadrangulable_geometry(GeometryData& gd);
	int df_check_ilp_side_subdivi(Eigen::VectorXi ilp_side_subdivi);

	//高山の残党
	template <typename PatchT>
	int df_quadrangulate(int simple_check, Eigen::VectorXi& l, PatchParam& param, PatchT& patch, std::vector<std::vector<Eigen::Vector2d>>& df_patches, int id, std::vector<std::vector<Eigen::VectorXi>> q_enum, int& now_pattern, int inc);
	template <typename PatchT>
	int df_quadrangulate_1(Eigen::VectorXi& l, PatchParam& param, PatchT& patch);
	template <typename PatchT>
	int df_quadrangulate_3(Eigen::VectorXi& l, PatchParam& param, PatchT& patch, std::vector<std::vector<Eigen::Vector2d>>& df_patches, int id, std::vector<std::vector<Eigen::VectorXi>> q_enum, int& now_pattern, int inc);
	template <typename PatchT>
	int df_quadrangulate_1_pentagon(Eigen::VectorXi& l, PatchParam& param, PatchT& patch);
	template <typename PatchT>
	int df_quadrangulate_3_6gon(Eigen::VectorXi& l, PatchParam& param, PatchT& patch, std::vector<std::vector<Eigen::VectorXi>> q_enum, int& now_pattern, int inc);    
}
