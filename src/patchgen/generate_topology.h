#include "decl.h"
#include <kt84/util.h>
#include <kt84/openmesh/flip_faces.hh>
#include "../patchgen_demo/curvenetwork/Patch.hh"
#include "../patchgen_demo/curvenetwork/decl.hh"

//dfrqd用
#include "ILP.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <kt84/eigen_def.hh>
#include <kt84/openmesh/edgeloop.hh>
#include <sstream>

template <typename PatchT>
int patchgen::df_quadrangulate(int simple_check, Eigen::VectorXi& l, PatchParam& param, PatchT& patch, std::vector<std::vector<Eigen::Vector2d>>& df_patches, int id, std::vector<std::vector<Eigen::VectorXi>> q_enum, int& now_pattern, int inc) {
	
	int num_sides = l.size();

	//① Simple Convex
	bool check1 = ((l.size() <= 5 && l.size() >= 3) && (l[0] != 0)) && (simple_check == 2);
	bool check3 = ((l.size() > 5) && (l[0] != 0)) && (simple_check == 2);
	if (check1) {
		patchgen::df_quadrangulate_1(l, param, patch);
		return 1;
	}
	else if(check3){
		patchgen::df_quadrangulate_3(l, param, patch, df_patches, id, q_enum, now_pattern, inc);//本当は4だったわ！
		return 3;
	}
	else {
		std::cout << "Type:No match\n";
	}

	return 0;
}

template <typename PatchT>
int patchgen::df_quadrangulate_1(Eigen::VectorXi& l, PatchParam& param, PatchT& patch) {
	std::cout << "Type:1 Simple Convex\n";

	if (l.size() == 3)
		std::cout << "Simple 3\n";
	else if (l.size() == 4)
		std::cout << "Simple 4\n";
	else if (l.size() == 5) 
		patchgen::df_quadrangulate_1_pentagon(l, param, patch);
	return 0;
}

template <typename PatchT>
int patchgen::df_quadrangulate_3(Eigen::VectorXi& l, PatchParam& param, PatchT& patch, std::vector<std::vector<Eigen::Vector2d>>& df_patches, int id, std::vector<std::vector<Eigen::VectorXi>> q_enum, int& now_pattern, int inc) {
	std::cout << "Type:3 Patches with |TV D|≦1\n";

	if (l.size() == 3)
		std::cout << "3-gon\n";
	else if (l.size() == 4)
		std::cout << "4-gon\n";
	else if (l.size() == 5)
		std::cout << "5-gon\n";
	else if (l.size() >= 6)
		df_quadrangulate_3_6gon(l, param, patch, q_enum, now_pattern, inc);
	else
		;
	return 0;
}

template <typename PatchT>
int patchgen::df_quadrangulate_1_pentagon(Eigen::VectorXi& l, PatchParam& param, PatchT& patch) {
	std::cout << "Simple Pentagon\n";
	
	int num_variables = l.size(); //ポリコードの数

	Eigen::MatrixXi A;
	Eigen::VectorXi X;
	A.resize(num_variables, num_variables);
	X.resize(num_variables, num_variables);
	//ポリコードは自分の今いる辺から前後に1つ飛ばしたところへそれぞれ通る
	/*A << 1, 1, 0, 0, 0,
		 0, 0, 1, 1, 0,
		 1, 0, 0, 0, 1,
		 0, 1, 0, 1, 0,
		 0, 0, 1, 0, 1;*/
	A/*^-1*/ << 1, 1, 1, -1, -1,
			   1, -1, -1, 1, 1,
			   1, 1, -1, -1, 1,
			   -1, 1, 1, 1, -1,
			   -1, -1, 1, 1, 1;
	X = A * l;
	X = X / 2;
	//std::cout << X << "\n";	
	/*
	patch.clear();
	typename PatchT::VHandle C[5];
	typename PatchT::VHandle E[6];
	for (int i = 0; i < 5; ++i) C[i] = add_tagged_vertex(patch, i, true);
	for (int i = 0; i < 6; ++i)	E[i] = add_tagged_vertex(patch, i, false);
	patch.add_face(E[5], E[0], C[0], E[1]);
	patch.add_face(E[5], E[1], C[1], E[2]);
	patch.add_face(E[5], E[2], C[2], E[3]);
	patch.add_face(E[5], E[3], C[3], E[4]);
	patch.add_face(E[5], E[4], C[4], E[0]);
	
	auto h_insert_a = patch.halfedge_handle(E[0]);
	if (X[4]-1 != 0)
		for (int i = 0; i < X[4]-1; ++i)
			insert_edgeloop(patch, h_insert_a);
	
	auto h_insert_b = patch.halfedge_handle(E[1]);
	if (X[1]-1 != 0)
		for (int i = 0; i < X[1]-1; ++i)
			insert_edgeloop(patch, h_insert_b);

	auto h_insert_c = patch.halfedge_handle(E[2]);
	if (X[2]-1 != 0)
		for (int i = 0; i < X[2]-1; ++i)
			insert_edgeloop(patch, h_insert_c);

	auto h_insert_d = patch.halfedge_handle(E[3]);
	if (X[0]-1 != 0)
		for (int i = 0; i < X[0]-1; ++i)
			insert_edgeloop(patch, h_insert_d);

	auto h_insert_e = patch.halfedge_handle(E[4]);
	if (X[3]-1 != 0)
		for (int i = 0; i < X[3]-1; ++i)
			insert_edgeloop(patch, h_insert_e);
	*/
	return 0;
}

template <typename PatchT>
int patchgen::df_quadrangulate_3_6gon(Eigen::VectorXi& l, PatchParam& param, PatchT& patch, std::vector<std::vector<Eigen::VectorXi>> q_enum, int& now_pattern, int inc) {
	int num_sides = l.size(); //ポリコードの数
	
	std::cout << num_sides << "-gon is input\n";
	std::cout << now_pattern << "\n";

	//6角形の場合，2つの五角形に分割されるとして，
	if (num_sides == 6)
	{
		int se = 1;//対象とするエッジ
		
		//n1 n2 n3 n4 x S1 S2 t1 = C
		int num_variables = 7;
		//std::vector<std::vector<Eigen::VectorXi>> q_enum;


		for (se = 0; se < num_sides; se++)
		//for (se = 0; se < num_sides / 2; se++)
		{
			std::vector<int> qe;
			q_enum.resize(num_sides);

			for (int n1 = 1; n1 < l[se % 6]; n1++)
			{
				for (int n3 = 1; n3 < l[(se + 3) % 6]; n3++)
				{
					for (int x = 1;; x++)
					{
						ILP ilp(num_variables);

						//サブパッチの合計が偶数
						int se1 = (se + 1) % 6;
						int se2 = (se + 2) % 6;
						int se3 = (se + 4) % 6;
						int se4 = (se + 5) % 6;
						ilp.add_constraint(kt84::make_Vector7d(1, 0, 1, 0, 1, -2, 0), EQ, -(l[se1] + l[se2]));
						ilp.add_constraint(kt84::make_Vector7d(0, 1, 0, 1, 1, 0, -2), EQ, -(l[se3] + l[se4]));

						//シンプルの条件
						ilp.add_constraint(kt84::make_Vector7d(-1, 0, 1, 0, 1, 0, 0), GE, -l[se1] + l[se2]);
						ilp.add_constraint(kt84::make_Vector7d(1, 0, 1, 0, 1, 0, 0), GE, l[se1] + l[se2]);
						ilp.add_constraint(kt84::make_Vector7d(1, 0, -1, 0, 1, 0, 0), GE, l[se1] - l[se2]);
						ilp.add_constraint(kt84::make_Vector7d(1, 0, -1, 0, -1, 0, 0), GE, -l[se1] - l[se2]);
						ilp.add_constraint(kt84::make_Vector7d(-1, 0, 1, 0, -1, 0, 0), GE, -l[se1] - l[se2]);

						ilp.add_constraint(kt84::make_Vector7d(0, -1, 0, 1, 1, 0, 0), GE, -l[se3] + l[se4]);
						ilp.add_constraint(kt84::make_Vector7d(0, 1, 0, 1, 1, 0, 0), GE, l[se3] + l[se4]);
						ilp.add_constraint(kt84::make_Vector7d(0, 1, 0, -1, 1, 0, 0), GE, l[se3] - l[se4]);
						ilp.add_constraint(kt84::make_Vector7d(0, 1, 0, -1, -1, 0, 0), GE, -l[se3] - l[se4]);
						ilp.add_constraint(kt84::make_Vector7d(0, -1, 0, 1, -1, 0, 0), GE, -l[se3] - l[se4]);

						//分割する場所の条件
						ilp.add_constraint(kt84::make_Vector7d(1, 1, 0, 0, 0, 0, 0), EQ, l[se % 6]);
						ilp.add_constraint(kt84::make_Vector7d(0, 0, 1, 1, 0, 0, 0), EQ, l[(se + 3) % 6]);

						//列挙のためのパラメータ調整
						ilp.add_constraint(kt84::make_Vector7d(1, 0, 0, 0, 0, 0, 0), GE, n1);
						ilp.add_constraint(kt84::make_Vector7d(0, 1, 0, 0, 0, 0, 0), GE, n1);
						ilp.add_constraint(kt84::make_Vector7d(0, 0, 1, 0, 0, 0, 0), GE, n3);
						ilp.add_constraint(kt84::make_Vector7d(0, 0, 0, 1, 0, 0, 0), GE, n3);
						ilp.add_constraint(kt84::make_Vector7d(0, 0, 0, 0, 1, 0, 0), GE, x);

						ilp.set_objective(kt84::make_Vector7d(0, 0, 0, 0, 1, 0, 0), false);

						if (!ilp.solve()) break;

						auto variables = ilp.get_variables();
						int tmp = 0;
						for (int i = 0; i < 5; ++i)
						{
							tmp += variables[i] * pow(10, (5 - i - 1) * 2);
						}
						qe.push_back(tmp);

						std::sort(qe.begin(), qe.end());
						qe.erase(std::unique(qe.begin(), qe.end()), qe.end());
					}
				}
			}
			q_enum[se].resize(qe.size());
			for (int n = 0; n < qe.size(); n++)
			{
				//std::cout << qe[n] << "\n";
				q_enum[se][n].resize(5);
				int tmp = qe[n];
				for (int m = 0; m < 5; m++)
				{
					int tmp2 = (int)(10 * (tmp / (int)pow(10, (5 - m) * 2)));
					tmp = (int)(tmp % (int)pow(10, (5 - m) * 2));
					tmp2 += (int)(tmp / (int)pow(10, (5 - m - 1) * 2));
					tmp = (int)(tmp % (int)pow(10, (5 - m - 1) * 2));
					q_enum[se][n][m] = tmp2;
					//std::cout << q_enum[se][n][m] << " ";
				}
				//std::cout << "\n";
			}
		}

		Eigen::Vector5i l1;
		Eigen::Vector5i l2;
		Eigen::MatrixXi A;
		Eigen::VectorXi X1;
		Eigen::VectorXi X2;
		A.resize(5, 5);
		X1.resize(5);
		X2.resize(5);

		for (;;){

			int count = 0;
			int select_i = 0;
			int select_j = 0;

			int max = 0;
			for (int i = 0; i < q_enum.size(); i++)
				max += q_enum[i].size();

			now_pattern = abs(now_pattern);
			now_pattern = now_pattern % max;

			for (int i = 0; i < (num_sides); i++)/////////回転によって計算のされ方が違う･･････＞＜
				//for (int i = 0; i < (num_sides / 2); i++)
			{
				for (int j = 0; j < q_enum[i].size(); j++)
				{
					for (int k = 0; k < 5; k++)
					{
						;//std::cout << q_enum[i][j][k] << " ";
					}
					if (now_pattern == count)
					{

						select_i = i;
						select_j = j;
						;//std::cout << "　←";
					}

					;//std::cout << "\n";
					count++;
				}
				;//std::cout << "\n";
			}
			
			//Eigen::Vector5i l1;
			//Eigen::Vector5i l2;
			
			int se1 = (select_i + 1) % 6;
			int se2 = (select_i + 2) % 6;
			int se3 = (select_i + 4) % 6;
			int se4 = (select_i + 5) % 6;
			l1 << q_enum[select_i][select_j][0], l[se1], l[se2], q_enum[select_i][select_j][2], q_enum[select_i][select_j][4];
			l2 << q_enum[select_i][select_j][3], l[se3], l[se4], q_enum[select_i][select_j][1], q_enum[select_i][select_j][4];
			A/*^-1*/ << 1, 1, 1, -1, -1,
				1, -1, -1, 1, 1,
				1, 1, -1, -1, 1,
				-1, 1, 1, 1, -1,
				-1, -1, 1, 1, 1;
			X1 = A * l1;
			X1 = X1 / 2;
			X2 = A * l2;
			X2 = X2 / 2;
			/*
			std::cout << "		";
			for (int p = 0; p < 5; p++)
				std::cout << X1[p] << " ";
			std::cout << "\n";
			std::cout << "		";
			for (int p = 0; p < 5; p++)
				std::cout << X2[p] << " ";
			std::cout << "\n";

			std::cout << "		l1 ";
			for (int p = 0; p < 5; p++)
				std::cout << l1[p] << " ";
			std::cout << "\n";
			std::cout << "		l2 ";
			for (int p = 0; p < 5; p++)
				std::cout << l2[p] << " ";
			std::cout << "\n";
			*/
			bool zero = (X1[0] > 0) && (X1[2] > 0) && (X1[3] != 0) && (X1[4] > 0)
				&& (X2[0] > 0) && (X2[2] > 0) && (X2[3] > 0) && (X2[4] > 0);
			bool edge_type = ((X1[1] + X2[1]) > 0);
			bool p0 = (X1[2] > X2[4]);
			bool p1 = (X1[2] == X2[4]);


			patch.clear();

			int n[6] = { 0, 1, 2, 3, 4, 5 };

			if (zero && edge_type)
			{
				if (p1)
				{
					std::cout << "pattern1" << "\n";//g.df_l << 3, 2, 2, 3, 2, 2;

					for (int i = 0; i < 6; i++)
						n[i] = ((n[i] - select_i) < 0) ? n[i] - select_i + 6 : n[i] - select_i;

					typename PatchT::VHandle C[6];
					typename PatchT::VHandle E[10];
					for (int i = 0; i < 6; ++i) C[i] = add_tagged_vertex(patch, i, true);
					for (int i = 0; i < 10; ++i) E[i] = add_tagged_vertex(patch, i, false);
					patch.add_face(E[0], E[8], E[7], C[n[0]]);
					patch.add_face(E[1], E[9], E[8], E[0]);
					patch.add_face(E[2], E[9], E[1], C[n[5]]);
					patch.add_face(E[3], E[9], E[2], C[n[4]]);
					patch.add_face(E[4], E[9], E[3], C[n[3]]);
					patch.add_face(E[5], E[8], E[9], E[4]);
					patch.add_face(E[6], E[8], E[5], C[n[2]]);
					patch.add_face(E[7], E[8], E[6], C[n[1]]);

					auto h_insert_a = patch.halfedge_handle(E[1]);
					if ((X1[1] + X2[1]) - 1 != 0)
						for (int i = 0; i < (X1[1] + X2[1]) - 1; ++i)
							insert_edgeloop(patch, h_insert_a);

					auto h_insert_b = patch.halfedge_handle(E[2]);
					if (X1[2] - 1 != 0)
						for (int i = 0; i < X1[2] - 1; ++i)
							insert_edgeloop(patch, h_insert_b);

					auto h_insert_c = patch.halfedge_handle(E[6]);
					if (X2[2] - 1 != 0)
						for (int i = 0; i < X2[2] - 1; ++i)
							insert_edgeloop(patch, h_insert_c);

					auto h_insert_d = patch.halfedge_handle(E[3]);
					if (X1[0] - 1 != 0)
						for (int i = 0; i < X1[0] - 1; ++i)
							insert_edgeloop(patch, h_insert_d);

					auto h_insert_e = patch.halfedge_handle(E[7]);
					if (X2[0] - 1 != 0)
						for (int i = 0; i < X2[0] - 1; ++i)
							insert_edgeloop(patch, h_insert_e);

					auto h_insert_f = patch.halfedge_handle(E[0]);
					if (X2[3] - 1 != 0)
						for (int i = 0; i < X2[3] - 1; ++i)
							insert_edgeloop(patch, h_insert_f);

					auto h_insert_g = patch.halfedge_handle(E[4]);
					if (X1[3] - 1 != 0)
						for (int i = 0; i < X1[3] - 1; ++i)
							insert_edgeloop(patch, h_insert_g);
				}
				else
				{
					std::cout << "pattern0" << "\n";
					typename PatchT::VHandle C[6];
					typename PatchT::VHandle E[14];
					for (int i = 0; i < 6; ++i) C[i] = add_tagged_vertex(patch, i, true);
					for (int i = 0; i < 14; ++i) E[i] = add_tagged_vertex(patch, i, false);


					for (int i = 0; i < 6; i++)
						n[i] = ((n[i] + select_i) > 5) ? n[i] + select_i - 6 : n[i] + select_i;


					if (p0)
					{
						std::cout << "true" << "\n";//g.df_l << 4, 3, 2, 4, 3, 2;

						patch.add_face(E[0], E[10], E[9], C[n[4]]);
						patch.add_face(E[1], E[13], E[10], E[0]);
						patch.add_face(E[2], E[13], E[1], C[n[5]]);
						patch.add_face(E[3], E[13], E[2], C[n[0]]);
						patch.add_face(E[4], E[11], E[13], E[3]);
						patch.add_face(E[5], E[11], E[4], C[n[1]]);
						patch.add_face(E[6], E[12], E[11], E[5]);
						patch.add_face(E[7], E[12], E[6], C[n[2]]);
						patch.add_face(E[8], E[12], E[7], C[n[3]]);
						patch.add_face(E[9], E[10], E[12], E[8]);
						patch.add_face(E[10], E[13], E[11], E[12]);

						auto h_insert_a = patch.halfedge_handle(E[9]);
						if ((X1[1] + X2[1]) - 1 != 0)
							for (int i = 0; i < (X1[1] + X2[1]) - 1; ++i)
								insert_edgeloop(patch, h_insert_a);

						auto h_insert_b = patch.halfedge_handle(E[5]);
						if (X2[4] - 1 != 0)
							for (int i = 0; i < X2[4] - 1; ++i)
								insert_edgeloop(patch, h_insert_b);

						auto h_insert_c = patch.halfedge_handle(E[0]);
						if (X1[4] - 1 != 0)
							for (int i = 0; i < X1[4] - 1; ++i)
								insert_edgeloop(patch, h_insert_c);

						auto h_insert_d = patch.halfedge_handle(E[3]);
						if (X2[3] - 1 != 0)
							for (int i = 0; i < X2[3] - 1; ++i)
								insert_edgeloop(patch, h_insert_d);

						auto h_insert_e = patch.halfedge_handle(E[1]);
						if ((X1[2] - X2[4]) - 1 != 0)
							for (int i = 0; i < (X1[2] - X2[4]) - 1; ++i)
								insert_edgeloop(patch, h_insert_e);

						auto h_insert_f = patch.halfedge_handle(E[2]);
						if (X2[0] - 1 != 0)
							for (int i = 0; i < X2[0] - 1; ++i)
								insert_edgeloop(patch, h_insert_f);

						auto h_insert_g = patch.halfedge_handle(E[8]);
						if (X1[3] - 1 != 0)
							for (int i = 0; i < X1[3] - 1; ++i)
								insert_edgeloop(patch, h_insert_g);

						auto h_insert_h = patch.halfedge_handle(E[7]);
						if (X1[0] - 1 != 0)
							for (int i = 0; i < X1[0] - 1; ++i)
								insert_edgeloop(patch, h_insert_h);
					}
					else
					{
						std::cout << "false" << "\n";//g.df_l << 4, 2, 3, 4, 2, 3;

						for (int i = 0; i < 6; i++)
							n[i] = (n[i] == 0) ? 5 : n[i] - 1;

						patch.add_face(E[0], E[10], E[9], C[n[4]]);
						patch.add_face(E[1], E[13], E[10], E[0]);
						patch.add_face(E[2], E[13], E[1], C[n[5]]);
						patch.add_face(E[3], E[13], E[2], C[n[0]]);
						patch.add_face(E[4], E[11], E[13], E[3]);
						patch.add_face(E[5], E[11], E[4], C[n[1]]);
						patch.add_face(E[6], E[12], E[11], E[5]);
						patch.add_face(E[7], E[12], E[6], C[n[2]]);
						patch.add_face(E[8], E[12], E[7], C[n[3]]);
						patch.add_face(E[9], E[10], E[12], E[8]);
						patch.add_face(E[10], E[13], E[11], E[12]);

						auto h_insert_a = patch.halfedge_handle(E[1]);
						if ((X1[1] + X2[1]) - 1 != 0)
							for (int i = 0; i < (X1[1] + X2[1]) - 1; ++i)
								insert_edgeloop(patch, h_insert_a);

						auto h_insert_e = patch.halfedge_handle(E[9]);
						if ((X2[4] - X1[2]) - 1 != 0)
							for (int i = 0; i < (X2[4] - X1[2]) - 1; ++i)
								insert_edgeloop(patch, h_insert_e);

						auto h_insert_c = patch.halfedge_handle(E[0]);
						if (X1[3] - 1 != 0)
							for (int i = 0; i < X1[3] - 1; ++i)
								insert_edgeloop(patch, h_insert_c);

						auto h_insert_b = patch.halfedge_handle(E[5]);
						if (X2[3] - 1 != 0)
							for (int i = 0; i < X2[3] - 1; ++i)
								insert_edgeloop(patch, h_insert_b);

						auto h_insert_d = patch.halfedge_handle(E[3]);
						if (X2[0] - 1 != 0)
							for (int i = 0; i < X2[0] - 1; ++i)
								insert_edgeloop(patch, h_insert_d);

						auto h_insert_g = patch.halfedge_handle(E[8]);
						if (X1[0] - 1 != 0)
							for (int i = 0; i < X1[0] - 1; ++i)
								insert_edgeloop(patch, h_insert_g);

						auto h_insert_f = patch.halfedge_handle(E[2]);
						if (X2[2] - 1 != 0)
							for (int i = 0; i < X2[2] - 1; ++i)
								insert_edgeloop(patch, h_insert_f);

						auto h_insert_h = patch.halfedge_handle(E[7]);
						if (X1[2] - 1 != 0)
							for (int i = 0; i < X1[2] - 1; ++i)
								insert_edgeloop(patch, h_insert_h);

					}

				}
				break;
			}
			else
			{
				if (inc == 0)
					now_pattern++;
				else if (inc == 1)
					now_pattern--;
				else
					now_pattern++;
			}
		}
		int c2 = 0;
		for (int i = 0; i < (num_sides); i++)/////////回転によって計算のされ方が違う･･････＞＜
			//for (int i = 0; i < (num_sides / 2); i++)
		{
			for (int j = 0; j < q_enum[i].size(); j++)
			{
				for (int k = 0; k < 5; k++)
				{
					std::cout << q_enum[i][j][k] << " ";
				}
				if (now_pattern == c2)
				{
					std::cout << "　←";
				}

				std::cout << "\n";
				c2++;
			}
			std::cout << "\n";
		}

		std::cout << "		";
		for (int p = 0; p < 5; p++)
			std::cout << X1[p] << " ";
		std::cout << "\n";
		std::cout << "		";
		for (int p = 0; p < 5; p++)
			std::cout << X2[p] << " ";
		std::cout << "\n";

	}
	return 0;
}

