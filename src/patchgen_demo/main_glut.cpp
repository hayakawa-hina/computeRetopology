#include <cstdlib>
#include <boost/lexical_cast.hpp>
#include <patchgen/decl.h>
#include <patchgen/generate_topology.h>
#include <kt84/graphics/graphics_util.hh>
#include <kt84/tw_util.h>
#include <kt84/glut_util.hh>
#include <kt84/util.h>
#include <kt84/eigen_util.hh>
#include <kt84/geometry/CameraFree.hh>
#include <kt84/MinSelector.hh>
#include "curvenetwork/Patch.hh"
#include "curvenetwork/decl.hh"

#include <kt84/openmesh/append_quad_strip.hh>
#include <kt84/openmesh/flip_faces.hh>


using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;


struct Globals {
    TwBar* bar = nullptr;
    CameraFree camera;
    
    curvenetwork::Patch patch;
	patchgen::GeometryData gd;

	//flag
	bool df_patch_flag = false;
	bool df_draw_flag = false;
	int df_draw_count = 0;
	bool df_quadrangulate_flag = false;
	bool df_solution_flag = false;

	//
	int num_sides = 5;
	int pattern_number = 0;//選択しているパターン(出力の)　⇒　GD
	int df_patch_id = 0;//選択しているパッチ(入力の)
	vector<vector<Vector2d>> df_patches;//最初のパッチ
	vector<Vector2d> df_draw_patch;
	vector<VectorXi> df_ls;
	Vector2d mouse_grobal;
    int selected_side = 0;
    int selected_variable = 0;
    curvenetwork::Patch::VHandle selected_vertex;
    int max_subdiv = 20;
	std::vector<std::vector<Eigen::VectorXi>> q_enum;
	Eigen::VectorXi change_subdivi;
	
	//std::vector<int> selected_patch;

    void init() {
        camera.eye << 0, 0, 5;
		
		int num_pattern = 8;

		df_patches.resize(num_pattern);

		df_patches[0].resize(5);
		df_patches[0][0] << 0.9, 0.3;//1.5→15
		df_patches[0][1] << 0.6, -1.0;//1.216→12
		df_patches[0][2] << -0.6, -1.0;//2.10→21
		df_patches[0][3] << -0.9, 0.3;//0.763→8
		df_patches[0][4] << 0, 1.2;
		
		df_patches[1].resize(6);
		df_patches[1][0] << 0.6, 0.5;//1.5→15
		df_patches[1][1] << 0.6, -1.0;//1.216→12
		df_patches[1][2] << -0.6, -1.2;//2.10→21
		df_patches[1][3] << -0.6, 0.9;//0.763→8
		df_patches[1][4] << -0.20, 1.55;//0.412→4
		df_patches[1][5] << 0.20, 1.45;//1.03→10
		
		df_patches[2].resize(6);
		df_patches[2][0] << 1.2, 0.0;//1.5→15
		df_patches[2][1] << 0.6, -1.42;//1.216→12
		df_patches[2][2] << -0.6, -1.42;//2.10→21
		df_patches[2][3] << -1.2, 0.0;//0.763→8
		df_patches[2][4] << -0.6, 1.42;//0.412→4
		df_patches[2][5] << 0.6, 1.42;//1.03→10

		df_patches[3].resize(5);
		df_patches[3][0] << 0.0, -1.5;//0.5→10
		df_patches[3][1] << -0.5, -1.5;//3.0→60
		df_patches[3][2] << -0.5, 1.5;//0.75→15
		df_patches[3][3] << 0.25, 1.5;//0.785→16
		df_patches[3][4] << 0.50, 0.755;//2.309→46
		
		df_patches[4].resize(7);
		df_patches[4][0] << 1.2, 0.0;//
		df_patches[4][1] << 0.6, -1.42;//
		df_patches[4][2] << -0.6, -1.42;//
		df_patches[4][3] << -1.2, -0.71;//
		df_patches[4][4] << -1.2, 0.71;//
		df_patches[4][5] << -0.6, 1.42;//
		df_patches[4][6] << 0.6, 1.42;//

		df_patches[5].resize(8);
		df_patches[5][0] << 1.2, -0.71;//
		df_patches[5][1] << 0.6, -1.42;//
		df_patches[5][2] << -0.6, -1.42;//
		df_patches[5][3] << -1.2, -0.71;//
		df_patches[5][4] << -1.2, 0.71;//
		df_patches[5][5] << -0.6, 1.42;//
		df_patches[5][6] << 0.6, 1.42;//
		df_patches[5][7] << 1.2, 0.71;//

		df_patches[6].resize(10);
		df_patches[6][0] << 1.00, 1.50;//0.781
		df_patches[6][1] << 1.60, 1.00;//1.011
		df_patches[6][2] << 1.75, 0.00;//1.75
		df_patches[6][3] << 1.75, -1.25;//0.269
		df_patches[6][4] << 1.65, -1.50;//3.8
		df_patches[6][5] << -1.65, -1.50;//0.269
		df_patches[6][6] << -1.75, -1.25;//1.75
		df_patches[6][7] << -1.75, 0.00;//1.011
		df_patches[6][8] << -1.60, 1.00;//0.781
		df_patches[6][9] << -1.00, 1.50;//2.5

		df_patches[7].resize(9);
		df_patches[7][0] << 0.75, 1.50;//0.781
		df_patches[7][1] << 1.25, 1.00;//1.011
		df_patches[7][2] << 1.75, 0.00;//1.75
		df_patches[7][3] << 1.00, -1.25;//0.269
		df_patches[7][4] << 0, -1.50;//3.8
		df_patches[7][5] << -1.00, -1.25;//1.75
		df_patches[7][6] << -1.75, 0.00;//1.011
		df_patches[7][7] << -1.25, 1.00;//0.781
		df_patches[7][8] << -0.75, 1.50;//2.5

		df_patches[7].resize(9);
		df_patches[7][0] << 0.75, 1.50;//0.781
		df_patches[7][1] << 1.25, 1.00;//1.011
		df_patches[7][2] << 1.75, 0.00;//1.75
		df_patches[7][3] << 1.00, -1.25;//0.269
		df_patches[7][4] << 0, -1.50;//3.8
		df_patches[7][5] << -1.00, -1.25;//1.75
		df_patches[7][6] << -1.75, 0.00;//1.011
		df_patches[7][7] << -1.25, 1.00;//0.781
		df_patches[7][8] << -0.75, 1.50;//2.5

		df_ls.resize(num_pattern);
		df_ls[0].resize(5);
		df_ls[0] << 2, 2, 3, 2, 5;
		df_ls[0] << 4, 4, 6, 4, 4;
		df_ls[1].resize(6);
		df_ls[1] << 15, 12, 22, 8, 4, 10;
		df_ls[2].resize(6);
		df_ls[2] << 2, 3, 3, 2, 3, 3;
		df_ls[2] << 4, 6, 2, 4, 6, 2;
		df_ls[3].resize(5);
		df_ls[3] << 10, 60, 15, 16, 48;
		//df_ls[3] << 4, 4, 4, 4, 4;
		//df_ls[3] << 10, 59, 15, 17, 49;
		//df_ls[3] << 20, 20, 20, 20, 20;
		df_ls[4].resize(7);
		df_ls[4] << 4, 2, 4, 3, 3, 4, 2;
		df_ls[4] << 3, 4, 2, 4, 2, 4, 3;
		df_ls[4] << 4, 4, 4, 4, 4, 4, 4;
		df_ls[4] << 3, 4, 2, 4, 2, 4, 3;
		df_ls[5].resize(8);
		df_ls[5] << 2, 4, 2, 4, 2, 4, 2, 4;
		df_ls[6].resize(10);
		df_ls[6] << 4, 4, 4, 4, 4, 4, 4, 4, 4, 4;
		df_ls[7].resize(9);
		df_ls[7] << 4, 4, 4, 4, 4, 4, 4, 4, 4;

		gd.corner_position = df_patches[0];
		gd.num_sides = df_ls[0].size();
		gd.side_subdivi = df_ls[0];
		gd.input_side_length = df_ls[0];

		change_subdivi.resize(gd.num_sides);
		change_subdivi = df_ls[0];
	}
} g;

void init_gl() {
    glClearColor(1, 1, 1, 1);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_POINT_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);    
    glHint(GL_POINT_SMOOTH_HINT, GL_NICEST);
    g.camera.auto_flip_y = false;
}

// >>AntTweakBar setup-------------------------------------------------
void init_bar() {
    g.bar = TwNewBar("patchgen");
    TwDefine("patchgen fontsize=3 size='300 320' valueswidth=fit color='153 204 51' alpha=200 text=dark");
    
	TwAddVarRW(g.bar, "dfRemesh_Patch", TW_TYPE_BOOLCPP, &g.df_patch_flag, "key=m");
	
	tw_util::AddButton(g.bar, "increment",
        [&](){
		if (g.df_patch_flag && !g.df_draw_flag){
			g.change_subdivi[g.selected_side] += 1;
			
			for (int i = 0; i < g.change_subdivi.size(); i++)
				std::cout << g.change_subdivi[i] << " ";
			std::cout << "\n";
		}
		
        }, "key=v");
    tw_util::AddButton(g.bar, "decrement",
        [&](){
		if (g.df_patch_flag && !g.df_draw_flag){
				g.change_subdivi[g.selected_side] -= 1;
				g.change_subdivi[g.selected_side] = (g.change_subdivi[g.selected_side] < 0) ? 0 : g.change_subdivi[g.selected_side];

				for (int i = 0; i < g.change_subdivi.size(); i++)
					std::cout << g.change_subdivi[i] << " ";
				std::cout << "\n";
			}
			
        }, "key=z");
    tw_util::AddButton(g.bar, "increment_by_2",
        [&](){
		if (g.df_patch_flag && !g.df_draw_flag){
				g.change_subdivi[g.selected_side] += 2;

				for (int i = 0; i < g.change_subdivi.size(); i++)
					std::cout << g.change_subdivi[i] << " ";
				std::cout << "\n";
			}
			
        }, "key=c");
    tw_util::AddButton(g.bar, "decrement_by_2",
        [&](){
		if (g.df_patch_flag && !g.df_draw_flag){
				g.change_subdivi[g.selected_side] -= 2;
				g.change_subdivi[g.selected_side] = (g.change_subdivi[g.selected_side] < 0) ? 0 : g.change_subdivi[g.selected_side];

				for (int i = 0; i < g.change_subdivi.size(); i++)
					std::cout << g.change_subdivi[i] << " ";
				std::cout << "\n";
			}
			
        }, "key=x");
	
    tw_util::AddButton(g.bar, "switch_pattern_forward",
        [&](){
		if (g.df_patch_flag && g.df_quadrangulate_flag && !g.df_draw_flag){
				g.patch.clear();
				
				g.gd.geometry_no++;
				if (g.gd.geometry_no >= g.gd.quadrangulable_geometry[g.gd.subdivi_edge_no].size()) {
					for (;;) {
						g.gd.subdivi_edge_no++;
						if (g.gd.subdivi_edge_no >= g.gd.quadrangulable_geometry.size()) {
							g.gd.subdivi_edge_no = 0;
						}
						if (g.gd.quadrangulable_geometry[g.gd.subdivi_edge_no].size() != 0) break;
					}
					g.gd.geometry_no = 0;
				}

				std::cout << "subdivi_edge_no = " << g.gd.subdivi_edge_no << "\n";
				std::cout << "geometry_no = " << g.gd.geometry_no << "\n";
			}
        }, "key=w");
    tw_util::AddButton(g.bar, "switch_pattern_backward",
        [&](){
		if (g.df_patch_flag && g.df_quadrangulate_flag && !g.df_draw_flag){
				g.patch.clear();

				g.gd.geometry_no--;
				if (g.gd.geometry_no < 0) {
					for (;;) {
						g.gd.subdivi_edge_no--;
						if (g.gd.subdivi_edge_no < 0) {
							g.gd.subdivi_edge_no = g.gd.quadrangulable_geometry.size() - 1;
						}
						if (g.gd.quadrangulable_geometry[g.gd.subdivi_edge_no].size() != 0) break;
					}
					g.gd.geometry_no = g.gd.quadrangulable_geometry[g.gd.subdivi_edge_no].size() - 1;
				}

				std::cout << "subdivi_edge_no = " << g.gd.subdivi_edge_no << "\n";
				std::cout << "geometry_no = " << g.gd.geometry_no << "\n";
			}
        }, "key=W");

	tw_util::AddButton(g.bar, "switch_Original_Patch",
		[&](){
		if (g.df_patch_flag && !g.df_draw_flag){
			g.patch.clear();

			g.df_patch_id = (g.df_patch_id == g.df_patches.size() - 1) ? 0 : g.df_patch_id + 1;
			
			g.pattern_number = 0;
			g.gd.input_side_length = g.df_ls[g.df_patch_id];
			g.gd.side_subdivi = g.gd.input_side_length;
			g.gd.corner_position = g.df_patches[g.df_patch_id];
			g.gd.num_sides = g.df_ls[g.df_patch_id].size();

			g.change_subdivi.resize(g.gd.num_sides);
			g.change_subdivi = g.df_ls[g.df_patch_id];
			g.df_quadrangulate_flag = false;
		}
	}, "key=n");

	tw_util::AddButton(g.bar, "Quadrangulate",
		[&](){
		if (g.df_patch_flag && !g.df_draw_flag){
			g.df_quadrangulate_flag = true;

			g.gd.ilp_side_subdivi = patchgen::df_calc_ilp_sides(g.change_subdivi);
			g.gd.side_subdivi = g.gd.ilp_side_subdivi;
			g.change_subdivi = g.gd.side_subdivi;

			patchgen::df_calc_quadrangulation(g.gd);
			g.df_solution_flag = (patchgen::df_check_ilp_side_subdivi(g.gd.ilp_side_subdivi) == 1) ? false : true;
			g.df_solution_flag = (patchgen::df_check_geometry(g.gd) == 1) ? false : g.df_solution_flag;

			if (g.df_solution_flag) {
				g.gd.geometry_no = 0;
				g.gd.subdivi_edge_no = 0;
				for (;;) {
					if (g.gd.quadrangulable_geometry[g.gd.subdivi_edge_no].size() != 0) break;
					g.gd.subdivi_edge_no++;
					if (g.gd.subdivi_edge_no >= g.gd.quadrangulable_geometry.size()) {
						g.gd.subdivi_edge_no = 0;
					}
				}
			}
		}
	}, "key=b");

	TwAddVarRW(g.bar, "Draw_Patch", TW_TYPE_BOOLCPP, &g.df_draw_flag, "key=m");

	tw_util::AddVarCB<int>(g.bar, "num_sides", TW_TYPE_INT32,
		[&](const int& value){
		g.num_sides = value;
	},
		[&](int& value){
		value = g.num_sides;
	}, "min=5 max=20 keyincr=UP keydecr=DOWN");
}
// <<AntTweakBar setup-------------------------------------------------

// >>GLUT callback functions---------------------------------------
void display_pre() {
    glut_util::defaultcb::display_pre();
    
    // set projection matrix
    double zNear = g.camera.eye.z() * 0.1;
    double zFar  = zNear * 100;
    double aspect_ratio = g.camera.width / static_cast<double>(g.camera.height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluPerspective(40, aspect_ratio, zNear, zFar);
    
    // set modelview matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
    gluLookAt(g.camera.eye, g.camera.center, g.camera.up);
}
void display_main() {
	auto& patch = g.patch;

	g.patch.clear();

	// draw patch パッチ内部の分割を描画
	if (g.df_patch_flag) {
		if (!g.df_draw_flag) {
			if (g.df_quadrangulate_flag) {
				if (g.df_solution_flag) {

					g.patch.df_set_geometry(g.gd);
					g.patch.df_draw();

				}
			}
			// draw boundary
			g.patch.df_draw_boundary(g.gd);
			// draw selected boundary
			g.patch.df_draw_selected_boundary(g.gd, g.selected_side);
			// draw number of edge subdivisions
			g.patch.df_draw_edge_subdivi(g.gd, g.change_subdivi);
			// draw singularities
			//g.patch.draw_singularities();
		}
		else {
				glLineWidth(4);
				::glColor3d(0, 0, 0);
				for (int i = 0; i < g.df_draw_patch.size(); i++) {
					glBegin(GL_LINE_STRIP);
					glVertex2d(g.df_draw_patch[i].x(), g.df_draw_patch[i].y());
					if (i == g.df_draw_patch.size() - 1)
						glVertex2d(g.mouse_grobal.x(), g.mouse_grobal.y());
					else
						glVertex2d(g.df_draw_patch[i + 1].x(), g.df_draw_patch[i + 1].y());
					glEnd();
				}
		}
	}
	
	// draw singularities
	if (g.df_patch_flag)
		;//patch.draw_singularities();//特異点の描画！！！！
}
void reshape(int width, int height) {
    glut_util::defaultcb::reshape(width, height);
    g.camera.reshape(width, height);
}
Vector2d get_mouse_pos(int x, int y) {
        Vector3d world_xyz = unproject(Vector3d(x, y, 0));
        Vector3d& camera_center = g.camera.center;
        Vector3d& camera_eye    = g.camera.eye;
        return camera_center.head(2) + (world_xyz - camera_center).head(2) * camera_eye.z() / (camera_eye.z() - world_xyz.z());
}
void mouse(int glut_button, int state, int x, int y) {
    if (TwEventMouseButtonGLUT(glut_button, state, x, y)) return glutPostRedisplay();
    y = g.camera.height - y;
    
    bool shift_pressed = (glutGetModifiers() & GLUT_ACTIVE_SHIFT) != 0;
    bool ctrl_pressed  = (glutGetModifiers() & GLUT_ACTIVE_CTRL ) != 0;
    bool alt_pressed   = (glutGetModifiers() & GLUT_ACTIVE_ALT  ) != 0;
    
    if (state == GLUT_UP) {
        if (g.camera.drag_mode != Camera::DragMode::NONE)
            return g.camera.mouse_up();
        else if (g.selected_vertex.is_valid())
            return g.selected_vertex.invalidate();
    }
    
    if (alt_pressed) {
        g.camera.mouse_down(x, y, ctrl_pressed ? Camera::DragMode::ZOOM : Camera::DragMode::PAN);
        return;
    }
    
    Vector2d mouse_pos = get_mouse_pos(x, y);

	g.mouse_grobal = mouse_pos;

	if (g.df_patch_flag){
		MinSelector<int> selected_side(-1);
		int size = g.df_patches[g.df_patch_id].size();
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < g.df_ls[g.df_patch_id][i]; ++j) {
				double t0 = i + j / static_cast<double>(g.df_ls[g.df_patch_id][i]);
				double t1 = i + (j + 1) / static_cast<double>(g.df_ls[g.df_patch_id][i]);
				Vector2d p0 = patchgen::df_get_boundary_geometry(size, t0, g.df_patches[g.df_patch_id]);
				Vector2d p1 = patchgen::df_get_boundary_geometry(size, t1, g.df_patches[g.df_patch_id]);
				auto dist = eigen_util::distance_to_line(p0, p1, mouse_pos, true);
				if (!dist) continue;
				selected_side.update(*dist, i);
			}
		}
		g.selected_side = selected_side.value;
	}

	if (g.df_draw_flag && g.df_patch_flag) {
		if (g.df_draw_count == 0) {
			g.df_draw_patch.clear();
		}
		
		if (state == GLUT_UP) {
			g.df_draw_patch.push_back(mouse_pos);
			g.df_draw_count++;
		}
		if (g.df_draw_count == g.num_sides) {
			g.df_draw_flag = false;
			g.df_draw_count = 0;
			
			g.df_patches.push_back(g.df_draw_patch);
			g.df_draw_patch.clear();
			VectorXi t(g.num_sides);
			for (int i = 0; i < g.num_sides; i++) {
				double norm = 0;
				if (i == g.num_sides - 1) {
					norm = (g.df_draw_patch[0] - g.df_draw_patch[g.num_sides - 1]).norm();
				}
				else {
					norm = (g.df_draw_patch[i] - g.df_draw_patch[i + 1]).norm();
				}
				t(i) = (int)(norm * 10);
			}
			g.df_ls.push_back(t);
		}
	}

    glutPostRedisplay();
}
void motion(int x, int y) {
    if (TwEventMouseMotionGLUT(x, y)) return glutPostRedisplay();
    y = g.camera.height - y;
    
    if (g.camera.drag_mode != Camera::DragMode::NONE) {
        g.camera.mouse_move(x, y);
        return glutPostRedisplay();
    }
    if (g.selected_vertex.is_valid()) {
        Vector2d mouse_pos = get_mouse_pos(x, y);
        g.patch.data(g.selected_vertex).laplaceDirect.value << mouse_pos, 0;
        return glutPostRedisplay();
    }
}
// <<GLUT callback functions---------------------------------------

int main_glut(int argc, char* argv[]) {
    glut_util::init(argc, argv, GLUT_DOUBLE | GLUT_RGBA, 0.8, true, "dfRemesh_demo",
        display_pre,
        display_main,
        glut_util::defaultcb::display_post,
        reshape,
        glut_util::defaultcb::keyboard,
        nullptr,
        glut_util::defaultcb::special,
        mouse,
        motion);
    
    // GLEW
    glewInit();
    
    // AntTweakBar
    TwInit(TW_OPENGL, NULL);
    TwGLUTModifiersFunc([](){return glutGetModifiers();});
    
    init_gl();
    init_bar();
    g.init();
    
    glutMainLoop();
    return 0;
}
