#include "decl.h"
#include <kt84/util.h>
#include <kt84/openmesh/flip_faces.hh>
#include "../patchgen_demo/curvenetwork/Patch.hh"
#include "../patchgen_demo/curvenetwork/decl.hh"
#include "ILP.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <kt84/eigen_def.hh>
#include <kt84/openmesh/edgeloop.hh>
#include <sstream>

Eigen::VectorXi patchgen::df_calc_ilp_sides(Eigen::VectorXi& l) {

	int num_sides = l.size();

	//ILP----------------------------
	std::cout << "----------ILP----------\n";
	//-------------------------------

	int num_variables = num_sides * 2 + 1;
	int C = 1;

	ILP ilp(num_variables);

	//�@��Βl�̐���
	Eigen::MatrixXd cons11;
	Eigen::MatrixXd cons12;
	Eigen::VectorXd cons11_rhs;
	Eigen::VectorXd cons12_rhs;
	cons11.resize(num_sides, num_variables);
	cons12.resize(num_sides, num_variables);
	cons11_rhs.resize(num_sides);
	cons12_rhs.resize(num_sides);
	for (int c = 0; c < num_sides; ++c) {
		for (int r = 0; r < num_sides; ++r) {
			if (r == c) {	//num_variables - ����slack �𔼕��ɕ����Čv�Z
				cons11(c, r) = 1.0 / l[c];
				cons11(c, r + num_sides) = -1;
				cons12(c, r) = 1.0 / l[c];
				cons12(c, r + num_sides) = 1;
			}
			else {
				cons11(c, r) = cons11(c, r + num_sides) = 0;
				cons12(c, r) = cons12(c, r + num_sides) = 0;
			}
		}
		cons11(c, num_variables - 1) = 0; //��������̂��߂�slack�ϐ�
		cons11_rhs(c) = 1; //�E�ӂ̒萔
		cons12(c, num_variables - 1) = 0;
		cons12_rhs(c) = 1;
	}

	//�A�􉽊w�I�����Ƃ̍��̐���
	//�B���������񕉂ł��鐧��
	//�ETVD��0�łȂ��Ƃ��C�e�T�C�h�̕�������2�ȏ�ł��鐧��
	Eigen::MatrixXd cons2;
	Eigen::VectorXd cons21_rhs;
	Eigen::VectorXd cons22_rhs;
	Eigen::VectorXd cons3_rhs;
	Eigen::VectorXd cons6_rhs;
	cons2.resize(num_sides, num_variables);
	cons21_rhs.resize(num_sides);
	cons22_rhs.resize(num_sides);
	cons3_rhs.resize(num_sides);
	cons6_rhs.resize(num_sides);
	for (int c = 0; c < num_sides; ++c) {
		for (int r = 0; r < num_sides; ++r) {
			if (r == c) {	//num_variables - ����slack �𔼕��ɕ����Čv�Z
				cons2(c, r) = 1;
				cons2(c, r + num_sides) = 0;
			}
			else {
				cons2(c, r) = cons2(c, r + num_sides) = 0;
			}
		}
		cons2(c, num_variables - 1) = 0; //��������̂��߂�slack�ϐ�
		cons21_rhs(c) = l[c] + C; //�E�ӂ̒萔
		cons22_rhs(c) = l[c] - C;
		cons3_rhs(c) = 1;
		cons6_rhs(c) = 2;
	}

	//�C�������̍��v�������ł��鐧��
	Eigen::VectorXd cons4;
	cons4.resize(num_variables);
	for (int r = 0; r < num_sides; ++r) {
		cons4(r) = 1;
		cons4(r + num_sides) = 0;
	}
	cons4(num_variables - 1) = -2;

	//�ړI�֐�
	Eigen::VectorXd object;
	object.resize(num_variables);
	for (int r = 0; r < num_sides; ++r) {
		object(r) = 0;
		object(r + num_sides) = 1;
	}
	object(num_variables - 1) = 0;

	//ILP�ɐ���@�`�C���Z�b�g
	ilp.add_constraint(cons11, LE, cons11_rhs);
	ilp.add_constraint(cons12, GE, cons12_rhs);
	ilp.add_constraint(cons2, LE, cons21_rhs);
	ilp.add_constraint(cons2, GE, cons22_rhs);
	ilp.add_constraint(cons2, GE, cons3_rhs);
	ilp.add_constraint(cons4, EQ, 0);

	//ILP�ɖړI�֐����Z�b�g
	ilp.set_objective(object, false);

	//ILP������(���ٓ_(v3-v5�y�A)�̐��͍l�����Ȃ�)
	if (!ilp.solve()) std::cout << "No solution (1)\n";
	auto variables1 = ilp.get_variables();

	//ILP�̖ړI�֐������Z�b�g
	ilp.refresh();

	//�D�T�C�h�̐��ɉ���������
	if (num_sides == 3) {
		int maxIndex;
		int max = l.maxCoeff(&maxIndex);
		Eigen::Vector7d cons51;
		for (int r = 0; r < num_sides; ++r) {
			cons51(r) = (r == maxIndex) ? -1 : 1;
			cons51(r + num_sides) = 0;
		}
		cons51(num_variables - 1) = 0;
		ilp.add_constraint(cons51, GE, 1);
	}
	else if (num_sides == 4) {
		Eigen::Vector9d cons52;
		Eigen::Vector9d cons53;
		cons52 << 1, 0, -1, 0, 0, 0, 0, 0, 0;
		cons53 << 0, 1, 0, -1, 0, 0, 0, 0, 0;
		ilp.add_constraint(cons52, EQ, 0);
		ilp.add_constraint(cons53, EQ, 0);
	}
	else if (num_sides == 5) {//���̃t�F�[�Y�ōŒ���2�ӂɑI�΂ꂽ�ӈȊO���ύX����邱�Ƃŋ��E��ɓ��ٓ_�����ꍇ����������I
		Eigen::Vector5d l_tmp;
		l_tmp << l[0] + l[1], l[1] + l[2], l[2] + l[3], l[3] + l[4], l[4] + l[0];
		int maxIndex;
		int max = l_tmp.maxCoeff(&maxIndex);

		Eigen::VectorXd cons54;
		cons54.resize(11);

		for (int r = 0; r < num_sides; ++r) {
			cons54(r + num_sides) = 0;
			if (r == maxIndex) {
				cons54(r) = -1;
				if (maxIndex == 4){
					cons54(0) = -1;
					cons54(r + num_sides) = 0;
				}
				else{
					++r;
					cons54(r) = -1;
					cons54(r + num_sides) = 0;
				}
			}
			else
				cons54(r) = 1;
		}
		cons54(num_variables - 1) = 0;

		ilp.add_constraint(cons54, GE, 1);
	}
	else if (num_sides >= 6) {//�T�u�p�b�`���V���v���ɂȂ��Ԃ��l�����ĕ����������肵�Ȃ���΂Ȃ�Ȃ��I
		Eigen::VectorXd l_tmp;
		l_tmp.resize(num_sides);
		for (int r = 0; r < num_sides; ++r)
			l_tmp(r) = (r != num_sides - 1) ? l[r] + l[r + 1] : l[r] + l[0];
		int maxIndex;
		int max = l_tmp.maxCoeff(&maxIndex);

		Eigen::VectorXd cons55;
		cons55.resize(num_variables);

		for (int r = 0; r < num_sides; ++r) {
			if (r == maxIndex || r == maxIndex + 1) {
				cons55(r) = -1;
				if (maxIndex == num_sides - 1)
					cons55(0) = -1;
			}
			else
				cons55(r) = 1;
			cons55(r + num_sides) = 0;
		}
		cons55(num_variables - 1) = 0;
		ilp.add_constraint(cons55, GE, 2 * abs(4 - num_sides));
	}

	//����E���Z�b�g
	ilp.add_constraint(cons2, GE, cons6_rhs);
	//ILP�ɖړI�֐����Z�b�g
	ilp.set_objective(object, false);

	//ILP������
	if (!ilp.solve()) std::cout << "No solution (2)\n";
	auto variables2 = ilp.get_variables();

	//�i�[����
	std::vector<Eigen::VectorXi> df_l;
	df_l.resize(3);
	df_l[0].resize(num_sides);
	df_l[1].resize(num_sides);
	df_l[2].resize(num_sides);

	for (int i = 0; i < num_sides; ++i) 
		df_l[0][i] = l[i];

	for (int i = 0; i < num_sides; ++i)
		df_l[1][i] = variables1(i);

	for (int i = 0; i < num_sides; ++i)
		df_l[2][i] = variables2(i);

	df_debug_calc_ilp(df_l);

	return df_l[2];
}

void patchgen::df_debug_calc_ilp(std::vector<Eigen::VectorXi> l)
{
	std::cout << "Geometry_length_input    = [";
	for (int i = 0; i < l[0].size(); ++i) {
		std::cout << " " << l[0][i];
	}
	std::cout << " ]\n";

	std::cout << "Subdivision_output   (1) = [";
	for (int i = 0; i < l[1].size(); ++i) {
		std::cout << " " << l[1][i];
	}
	std::cout << " ]\n";

	std::cout << "                     (2) = [";
	for (int i = 0; i < l[2].size(); ++i) {
		std::cout << " " << l[2][i];
	}
	std::cout << " ]\n";
}