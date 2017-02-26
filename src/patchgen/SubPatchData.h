#pragma once
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <vector>
#include <string>
#include <utility>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

namespace patchgen {

	struct IdxCounter;
	struct SubPatchEdgeData;
	struct SubPatchData;
	struct OverHexagonILP;

	struct IdxCounter {
		//SubPatchDataçÏê¨ÇÃÇΩÇﬂÇÃidx
		int SubPatchEdgeData_idx;
		int SubPatchData_idx;

		IdxCounter() {
			SubPatchEdgeData_idx = SubPatchData_idx = 0;
		}

		int getPatchIdx() {	return SubPatchData_idx++; }
		int getEdgeIdx() { return SubPatchEdgeData_idx++; }
	};

	struct OverHexagonILP {
		int num_variables;
		std::vector<Eigen::Vector2i> ilp_counter;
		std::vector<Eigen::Vector2i> ilp_counter_slack;

		Eigen::MatrixXd cons4;
		Eigen::VectorXd cons4_rhs;
		Eigen::MatrixXd cons0;
		Eigen::VectorXd cons0_rhs;
		Eigen::MatrixXd cons5;
		Eigen::VectorXd cons5_rhs;
		Eigen::MatrixXd cons3;
		Eigen::VectorXd cons3_rhs;
		Eigen::VectorXd cons0_ns_rhs;

		std::vector<SubPatchData> s_leaf;

		int first_subdivi_edge;
		std::vector<int> selected_subdivi_edge;

		int ilp_loop_count;
		int ilp_solve_count;
	};

	struct SubPatchEdgeData {
		int idx;
		int num_subdivi;
		int ilp_variable_idx;
		SubPatchData *belong_patch;

		SubPatchEdgeData *p_pointer;
		SubPatchEdgeData *c_pointer_r;
		SubPatchEdgeData *c_pointer_l;
		SubPatchEdgeData *same_edge_pointer;
		OpenMesh::PolyConnectivity::EHandle e_hande;
		OpenMesh::PolyConnectivity::VHandle v_handle_from;
		OpenMesh::PolyConnectivity::VHandle v_handle_to;


		//std::vector<int> ilp_num_subdivi;
		int ilp_num_subdivi;
		OpenMesh::PolyConnectivity::HHandle polyline;
		int num_polychord;

		SubPatchEdgeData() {
			ilp_variable_idx = -1;
		}

		~SubPatchEdgeData() {
			
		}

		void deleteSubPatchData() {
			
			std::cout << "delete_edge_id = " << idx << "\n";
			
		}
	};

	struct SubPatchData {
		int idx;
		int num_sides;
		std::vector<SubPatchEdgeData*> side_subdivi;

		SubPatchData *p_pointer;
		SubPatchData *c_pointer_r;
		SubPatchData *c_pointer_l;
		OpenMesh::PolyConnectivity::FHandle sub_face;

		SubPatchData() {
			;
		}

		~SubPatchData () {
			;
		}

		void deleteSubPatchData() {
			
			std::cout << "delete_patch_id = " << idx << "\n";
			for (int i = 0; i < side_subdivi.size(); i++) {
				delete side_subdivi[i];
			}			
		}
	};
}