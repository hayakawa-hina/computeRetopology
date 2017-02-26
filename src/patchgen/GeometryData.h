#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <utility>
#include "SubPatchData.h"

namespace patchgen {

	struct GeometryData {

		int num_sides;
		bool is_simple;
		
		std::vector<Eigen::Vector2d> corner_position;

		Eigen::VectorXi input_side_length;
		//std::vector<Eigen::VectorXi> ilp_side_subdivi;
		Eigen::VectorXi ilp_side_subdivi;


		Eigen::VectorXi side_subdivi;

		//�p�b�`����������Ƃ��̏��
		Eigen::VectorXi patch_geometry;
		std::vector<SubPatchData> subpatch;
		std::vector<std::vector<SubPatchData>> subpatch_leaf;
		std::vector<std::vector<Eigen::VectorXi>> quadrangulable_geometry;

		//�I�����Ă���p�b�`
		int subdivi_edge_no;
		int geometry_no;

		GeometryData() {
			;
		}
	};
}