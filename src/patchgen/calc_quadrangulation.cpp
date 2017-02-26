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


void patchgen::df_calc_quadrangulation(GeometryData& gd)
{
	std::cout << "----------Quadrangulation_start----------\n";

	if (df_check_quadrangulation_type(gd) == 0)
		std::cout << "----------Quadrangulation_end----------\n";
	else
		std::cout << "----------Quadrangulation_error----------\n";
}

int patchgen::df_check_quadrangulation_type(GeometryData& gd)
{
	int num_sides = gd.num_sides;

	if (num_sides == 3) return 1;
	else if (num_sides == 4) return 1;
	else if (num_sides == 5) return patchgen::df_set_simplepatch(gd);
	else if (num_sides >= 6) return patchgen::df_over_hexagon(gd);
	else if (num_sides <= 2) return 1;

	return -1;
}