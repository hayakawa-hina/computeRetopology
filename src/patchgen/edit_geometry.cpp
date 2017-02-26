#include "decl.h"
using namespace std;
using namespace Eigen;

Vector2d patchgen::df_get_boundary_geometry(int num_sides, long double t, vector<Vector2d> df_patch)
{

	if (t == num_sides) return df_patch[0];

	int i1 = static_cast<int>(t);
	int i2 = (i1 + 1) % num_sides;
	long double s = t - i1;
	return (1 - s) * df_patch[i1] + s * df_patch[i2];
}

