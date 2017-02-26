#include "Vertex.hh"
#include "Halfedge.hh"
#include "Halfchain.hh"
#include "Edgechain.hh"
#include "Patch.hh"
#include "Core.hh"
#include "helper.hh"
#include "Circulator.hh"
#include <kt84/util.h>
#include <kt84/eigen_util.hh>
#include <kt84/openmesh/append_quad_strip.hh>
#include <kt84/openmesh/flip_faces.hh>
#include <kt84/graphics/phong_tessellation.hh>
#include <kt84/graphics/graphics_util.hh>
#include <patchgen/decl.h>
#include <patchgen/generate_topology.h>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/PolyConnectivity.hh>
#include <OpenMesh/Core/Mesh/Traits.hh>
#include <OpenMesh/Core/Mesh/FinalMeshItemsT.hh>
#include <OpenMesh/Core/Mesh/AttribKernelT.hh>
#include <OpenMesh/Core/Mesh/PolyMeshT.hh>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>
#include <boost/lexical_cast.hpp>
#include <kt84/glut_util.hh>
using namespace std;
using namespace Eigen;
using namespace kt84;
using namespace kt84::graphics_util;

curvenetwork::Patch::Patch(int id_)
    : id        (id_)
    , is_deleted()
    , flag      ()
{}

void curvenetwork::Patch::clear() {
    PatchBase::clear();
    invalidate_displist();
}
void curvenetwork::Patch::set_halfchain(Halfchain* halfchain_) {
    halfchain = halfchain_;
    for (auto c = halfchain; ; ) {
        c->patch = this;
        c = c->next();
        if (c == halfchain)
            return;
    }
}
curvenetwork::Patch::HHandle curvenetwork::Patch::opposite_boundary_halfedge(HHandle boundary_halfedge) const {
    auto& h_data = data(boundary_halfedge);
    auto c = h_data.halfchain;
    auto c_opposite = c->opposite();
    
    auto p_opposite = opposite_patch(boundary_halfedge);
    if (!p_opposite)
        return HHandle();
    
    for (auto h : p_opposite->halfedges()) {
        auto& h_opposite_data = p_opposite->data(h);
        if (h_opposite_data.halfchain == c_opposite && h_opposite_data.index_wrt_halfchain == c->num_subdiv() - 1 - h_data.index_wrt_halfchain)
            return h;
    }
    
    return HHandle();
}
curvenetwork::Patch*         curvenetwork::Patch::opposite_patch            (HHandle boundary_halfedge) const {
    auto& h_data = data(boundary_halfedge);
    auto c = h_data.halfchain;
    auto c_opposite = c->opposite();
    return c_opposite->patch;
}
curvenetwork::Vertex*        curvenetwork::Patch::vertex_patch_to_curveNetwork(VHandle v) const {
    if (!is_boundary(v))
        return nullptr;
    
    auto h = prev_halfedge_handle(halfedge_handle(v));
    if (data(h).index_wrt_halfchain != 0)
        return nullptr;
    
    return data(h).halfchain->halfedge_front->from_vertex();
}
curvenetwork::Patch::VHandle curvenetwork::Patch::vertex_curveNetwork_to_patch(curvenetwork::Vertex* v) const {
    if (!v->is_corner())
        return VHandle();
    
    for (auto patch_v : vertices()) {
        if (!is_boundary(patch_v))
            continue;
        
        auto w = vertex_patch_to_curveNetwork(patch_v);
        if (w == v)
            return patch_v;
    }
    return VHandle();
}
curvenetwork::Patch::VHandle curvenetwork::Patch::patch_corner(int corner_index) const {
    for (auto v : vertices()) {
        if (data(v).patchgen.corner_index == corner_index) return v;
    }
    return VHandle();
}
pair<curvenetwork::Patch::HHandle, curvenetwork::Patch::HHandle> curvenetwork::Patch::get_side_halfedges(int side_index) const {
    auto h_front = prev_halfedge_handle(halfedge_handle(patch_corner(side_index)));
    auto h = h_front;
    while (data(from_vertex_handle(h)).patchgen.corner_index == -1)
        h = prev_halfedge_handle(h);
    auto h_back = h;
    return make_pair(h_front, h_back);
}
Polyline_PointNormal  curvenetwork::Patch::get_side_polyline (int side_index) const {
    Polyline_PointNormal polyline;
    auto h_pair = get_side_halfedges(side_index);
    auto h = h_pair.first;
    polyline.push_back(data(to_vertex_handle(h)).pn);
    while (true) {
        polyline.push_back(data(from_vertex_handle(h)).pn);
        if (h == h_pair.second) break;
        h = prev_halfedge_handle(h);
    }
    return polyline;
}
double curvenetwork::Patch::distance(const Vector3d& p) const {
    double dist_min = util::dbl_max();
    for (auto f : faces()) {
        auto fv = cfv_iter(f);
        Vector3d p0 = data(*fv).pn.head(3);    ++fv;
        Vector3d p1 = data(*fv).pn.head(3);    ++fv;
        Vector3d p2 = data(*fv).pn.head(3);    ++fv;
        Vector3d p3 = data(*fv).pn.head(3);
        double dist0 = *eigen_util::distance_to_triangle(p0, p1, p2, p, true);
        double dist1 = *eigen_util::distance_to_triangle(p0, p2, p3, p, true);
        dist_min = util::min(dist_min, dist0, dist1);
    }
    return dist_min;
}


void curvenetwork::Patch::add_padding(const patchgen::PatchParam& param) {
    const int num_sides = param.get_num_sides();
    // IMPORTANT NOTE: boundary halfedges are oriented CLOCKWISE!
    Patch::HHandle h_boundary_front;
    for (auto v : vertices()) {
        if (data(v).patchgen.corner_index == 0) {
            h_boundary_front = halfedge_handle(v);
            break;
        }
    }
    
    for (int i = num_sides - 1; i >= 0; --i) {
        Patch::HHandle h_boundary_back = h_boundary_front;
        while (data(to_vertex_handle(h_boundary_back)).patchgen.corner_index == -1)
            h_boundary_back = next_halfedge_handle(h_boundary_back);
        
        for (int j = 0; j < param.p[i]; ++j) {
            // clear corner flag for current corner vertices
            data(from_vertex_handle(h_boundary_front)).patchgen.corner_index = -1;
            data(to_vertex_handle  (h_boundary_back )).patchgen.corner_index = -1;
            
            kt84::append_quad_strip(*this, h_boundary_front, h_boundary_back);
            
            // set corner flag for new corner vertices
            data(from_vertex_handle(h_boundary_front)).patchgen.corner_index = (i + 1) % num_sides;
            data(to_vertex_handle  (h_boundary_back )).patchgen.corner_index = i;
        }
        
        // go to next (clockwise) side
        h_boundary_front = next_halfedge_handle(h_boundary_back);
    }
    assert(data(from_vertex_handle(h_boundary_front)).patchgen.corner_index == 0);
}

/*
patchgen::PatchParam curvenetwork::Patch::generate_topology(const Eigen::VectorXi& num_subdiv) {
    // get default parameter
    auto param = patchgen::get_default_parameter(num_subdiv);
    // generate topology
    generate_topology(param);
	return param;
}
*/
/*
void curvenetwork::Patch::generate_topology(const patchgen::PatchParam& param) {
    //patchgen::generate_topology(param, *this);
    garbage_collection();
#ifndef NDEBUG
    debugInfo_get(false, false, false);
#endif
    laplaceSolver_init();
    set_halfedgeData();
    set_boundary_pn();
}
*/

void curvenetwork::Patch::set_halfedgeData() {
	std::cout << "ok";
    auto v_corner = patch_corner(0);
    if (!v_corner.is_valid())
        // bug!
        return;
    auto h = prev_halfedge_handle(halfedge_handle(v_corner));
	std::cout << "ok";
	/*
	auto c = halfchain;
    if (!c)
        // bug!
        return;	
	std::cout << "ok";
    */
	
	// debug!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//#ifndef NDEBUG
    cout << "hogehoge!=============\n";
    for (auto h2 = h; ; ) {
        cout << "halfedge: " << h2.idx() << ", to-vertex corner: " << data(to_vertex_handle(h2)).patchgen.corner_index << endl;
        h2 = prev_halfedge_handle(h2);
        if (h2 == h)
            break;
    }
	/*
    for (auto c2 = c; ; ) {
        cout << "halfchain: " << c2->id << ", num_subdiv: " << c2->num_subdiv() << endl;
        c2 = c2->next();
        if (c2 == c)
            break;
    }
	*/
    cout << "\n\n";
//#endif
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
	/*
    while (true) {
        int index_wrt_halfchain = 0;
        
        for (int i = 0; i < c->num_subdiv(); ++i) {
            data(h).halfchain = c;
            data(h).index_wrt_halfchain = index_wrt_halfchain;
            
            h = prev_halfedge_handle(h);
            ++index_wrt_halfchain;
        }
        
        c = c->next();
        if (c == halfchain) {
            assert(data(to_vertex_handle(h)).patchgen.corner_index == 0);
            break;
        }
    }
	*/
}



void curvenetwork::Patch::laplaceSolver_init() {
    for (auto v : vertices()) {
        auto& vdata = data(v);
        
        vdata.laplaceDirect.is_fixed = vdata.laplaceIterative.is_fixed = is_boundary(v);
    }

    laplaceDirect_factorize();
}


void curvenetwork::Patch::set_boundary_pn() {
    auto v_corner = patch_corner(0);
    if (!v_corner.is_valid())
        // bug!
        return;
    auto h0 = prev_halfedge_handle(halfedge_handle(v_corner));
    for (auto h = h0; ; ) {
        auto& hdata = data(h);
        auto& pn = data(to_vertex_handle(h)).pn;
        auto c = hdata.halfchain;
        auto h_opposite = opposite_boundary_halfedge(h);
        Patch* opposite_patch = c->opposite()->patch;
        if (h_opposite.is_valid() && opposite_patch->flag == 0) {       // opposite_patch->flag is 1 if its vertex position is invalid
            // copy from adjacent patch boundary
            pn = opposite_patch->data(opposite_patch->from_vertex_handle(h_opposite)).pn;
        } else {
            // sample pn on Halfchain using arc-length parameterization
            double t = hdata.index_wrt_halfchain / static_cast<double>(c->num_subdiv());
            pn = c->pn_at(t);
        }
        h = prev_halfedge_handle(h);
        if (h == h0) break;
    }
    flag = 0;       // mark this patch as having valid vertex positions
}

void curvenetwork::Patch::render_phong_fill() {
    displist_phong_fill.render([&] () {
        for (auto f : faces()) {
            phong_tessellation::begin(phong_tessellation::Mode::TRIANGLE_FAN);
            for (auto v : fv_range(f))
                phong_tessellation::vertex(data(v).pn);
            phong_tessellation::end();
        }
    });
}
void curvenetwork::Patch::render_phong_line_interior() {
    displist_phong_line_interior.render([&] () {
        // render all but boundary edges
        phong_tessellation::begin(phong_tessellation::Mode::LINES);
        for (auto e : edges()) {
            if (is_boundary(e))
                continue;
            
            auto v0 = to_vertex_handle(halfedge_handle(e, 0));
            auto v1 = to_vertex_handle(halfedge_handle(e, 1));
            phong_tessellation::vertex(data(v0).pn);
            phong_tessellation::vertex(data(v1).pn);
        }
        phong_tessellation::end();
    });
}
void curvenetwork::Patch::render_phong_line_boundary() {
    displist_phong_line_boundary.render([&] () {
        // boundary edges only
        phong_tessellation::begin(phong_tessellation::Mode::LINES);
        for (auto e : edges()) {
            if (!is_boundary(e))
                continue;
            
            auto v0 = to_vertex_handle(halfedge_handle(e, 0));
            auto v1 = to_vertex_handle(halfedge_handle(e, 1));
            phong_tessellation::vertex(data(v0).pn);
            phong_tessellation::vertex(data(v1).pn);
        }
        phong_tessellation::end();
    });
}
void curvenetwork::Patch::render_flat_fill() {
    displist_flat_fill.render([&] () {
        glBegin(GL_QUADS);
        for (auto f : faces()) {
            Vector3d normal = Vector3d::Zero();
            for (auto v : fv_range(f))
                normal += data(v).pn.tail(3);
            normal.normalize();
            glNormal3d(normal);
            for (auto v : fv_range(f))
                glVertex3dv(&data(v).pn[0]);
        }
        glEnd();
    });
}
void curvenetwork::Patch::render_flat_line_interior() {
    displist_flat_line_interior.render([&] () {
        glBegin(GL_LINES);
        for (auto e : edges()) {
            if (is_boundary(e))
                continue;
            for (int i = 0; i < 2; ++i)
                glVertex3dv(&data(to_vertex_handle(halfedge_handle(e, i))).pn[0]);
        }
        glEnd();
    });
}
void curvenetwork::Patch::render_flat_line_boundary() {
    displist_flat_line_boundary.render([&] () {
        glBegin(GL_LINES);
        for (auto e : edges()) {
            if (!is_boundary(e))
                continue;
            for (int i = 0; i < 2; ++i)
                glVertex3dv(&data(to_vertex_handle(halfedge_handle(e, i))).pn[0]);
        }
        glEnd();
    });
}
void curvenetwork::Patch::render_irregular_vertices() const {
    static Vector3d color_table[10] = {
        { 0.0, 0.0, 0.0 },
        { 0.0, 0.0, 0.0 },
        { 0.1, 0.3, 1.0 },
        { 0.2, 0.6, 1.0 },
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.6, 0.0 },
        { 1.0, 0.3, 0.0 },
        { 0.7, 0.0, 0.9 },
        { 0.4, 0.8, 0.3 },
        { 0.0, 0.5, 0.0 },
    };
    glBegin(GL_POINTS);
    for (auto v : vertices()) {
        int valence_ = valence(v);
        bool is_corner = data(v).patchgen.corner_index != -1;
        if (is_corner)
            valence_ += 2;
        else if (is_boundary(v))
            ++valence_;
        if (valence_ == 4)
            continue;
        glColor3d(color_table[valence_]);
        glVertex3dv(&data(v).pn[0]);
    }
    glEnd();
}
void curvenetwork::Patch::invalidate_displist() {
    displist_phong_fill         .invalidate();
    displist_phong_line_interior.invalidate();
    displist_phong_line_boundary.invalidate();
    displist_flat_fill          .invalidate();
    displist_flat_line_interior .invalidate();
    displist_flat_line_boundary .invalidate();
}
/*
void curvenetwork::Patch::draw_singularities()
{
	for (auto v : vertices()) {
		int valence = 0;
		for (auto vv = cvv_iter(v); vv.is_valid(); ++vv, ++valence);
		if (is_boundary(v)) ++valence;
		if (data(v).patchgen.corner_index != -1) ++valence;
		if (valence == 4) continue;
		assert(valence == 3 || valence == 5);
		Vector3d p = data(v).laplaceDirect.value;
		Vector3d color = valence == 3 ? Vector3d(0, 0.6, 0.9) : Vector3d(1, 0.7, 0);
		glColor3d(0, 0, 0);    glPointSize(12);   glBegin(GL_POINTS);    glVertex3d(p);    glEnd();
		glColor3d(color);      glPointSize(10);   glBegin(GL_POINTS);    glVertex3d(p);    glEnd();
	}
}
*/
void curvenetwork::Patch::df_set_geometry(patchgen::GeometryData gd) 
{
	int num_sides = gd.corner_position.size();

	df_set_base_patch(gd);
	


	curvenetwork::Patch::HHandle h;
	curvenetwork::Patch::VHandle v0;
	v0 = gd.subpatch[gd.subdivi_edge_no].side_subdivi[0]->v_handle_from;
	for (auto f : faces()) {
		for (auto vv = cfh_iter(f); vv.is_valid(); ++vv) {
			if (from_vertex_handle(vv.handle()).idx() == v0.idx()) {
				h = vv.handle();
				break;
			}
		}
	}

	for (int i = 0; i < num_sides; ++i) {
		for (int j = 0; j < gd.side_subdivi[i]; ++j) {
			auto& vdata = data(from_vertex_handle(h)).laplaceDirect;
			long double t = i + j / static_cast<long double>(gd.side_subdivi[i]);
			vdata.value << patchgen::df_get_boundary_geometry(num_sides, t, gd.corner_position), 0;
			vdata.is_fixed = true;
			h = next_halfedge_handle(h);
		}
	}

	laplaceDirect_factorize();
	laplaceDirect_solve();


}

void curvenetwork::Patch::df_draw() 
{
	//FACE
	::glColor3d(0.7, 0.7, 0.7);
	glBegin(GL_POLYGON);
	for (auto f : faces()){
		for (auto v = cfv_iter(f); v.is_valid(); ++v){
			glVertex3d(data(*v).laplaceDirect.value);
		}
	}
	glEnd();

	//EDGE
	for (auto v : vertices()){
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glLineWidth(2);
		::glColor3d(0, 0, 0);
		glBegin(GL_LINE_STRIP);
		for (auto vv = cvv_iter(v); vv.is_valid(); ++vv)
		{
			glVertex3d(data(*vv).laplaceDirect.value);
			glVertex3d(data(v).laplaceDirect.value);
		}
		glEnd();
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}


void curvenetwork::Patch::df_debug_vertices_data()
{
	int i = 0;
	for (auto v : vertices())
	{
		i++;
		std::cout << "vertex_index = " << v.idx() << " : position = " << data(v).laplaceDirect.value.x() << " " << data(v).laplaceDirect.value.y() << "\n";

		for (auto vv = cvv_iter(v); vv.is_valid(); ++vv)
			std::cout << "		vertex_index = " << vv.handle().idx() << "\n";

		std::cout << "			halfedge_handle = " << halfedge_handle(v) << "\n";

		for (auto vv = cve_iter(v); vv.is_valid(); ++vv)
			std::cout << "		edge_index = " << vv.handle().idx() << "\n";
	}
	std::cout << "num_vertices " << i << "\n";
}

void curvenetwork::Patch::df_debug_faces_data()
{
	int j = 0, k = 0;
	for (auto f : faces()) {
		j++;
		k = 0;
		for (auto vv = cfh_iter(f); vv.is_valid(); ++vv) {
			std::cout << "		halfedge_index = " << vv.handle().idx() << "\n";
			k++;
		}
		std::cout << "	num_faces_halfedges = " << k << "\n";
	}
	std::cout << "num_faces " << j << "\n";
}

void curvenetwork::Patch::df_debug_edges_data()
{
	int k = 0;
	for (auto e : edges())
	{
		k++;
	}
	std::cout << "num_edges " << k << "\n";
}

void curvenetwork::Patch::df_debug_halfedges_data()
{
	int l = 0;
	for (auto h : halfedges()){
		l++;

		std::cout << "halfedge_index = " << h.idx() << "\n";
		std::cout << "	next_halfedge_index = " << next_halfedge_handle(h).idx() << " ";
		std::cout << "  prev_halfedge_index = " << prev_halfedge_handle(h).idx() << "\n";
		std::cout << "		from_vertex_index = " << from_vertex_handle(h).idx() << " ";
		std::cout << "		to  _vertex_index = " << to_vertex_handle(h).idx() << "\n";

	}
	std::cout << "num_halfedges " << l << "\n";
}


std::vector<curvenetwork::Patch::HHandle> curvenetwork::Patch::df_find_polyline_halfedge(curvenetwork::Patch::HHandle h)
{
	std::vector<curvenetwork::Patch::HHandle> fh;
	fh.resize(2);

	fh[0] = h;
	fh[1] = opposite_halfedge_handle(h);

	for (;;)
	{
		if (!is_boundary(fh[0])) break;
		fh[0] = opposite_halfedge_handle(next_halfedge_handle(next_halfedge_handle(fh[0])));
	}
	for (;;)
	{
		if (!is_boundary(fh[1])) break;
		fh[1] = opposite_halfedge_handle(next_halfedge_handle(next_halfedge_handle(fh[1])));
	}

	//std::cout << h.idx() << " -> " << fh[0].idx() << ", " << fh[1].idx() << "\n";

	return fh;
}


void curvenetwork::Patch::df_set_base_patch(patchgen::GeometryData gd)
{
	int num_sides = gd.num_sides;
	//-----------特異点が2つ以上のパターン（シンプルでない，六角形以上）はパッチの分割を敢行
	//-----------エッジに特異点が乗っている場合は基本形の頂点数を変更！

	if (false)
	{
		;
	}
	else if (num_sides >= 5) {
		int num_v_base = num_sides;

		df_set_subpatch_over_hexagon(&gd.subpatch[gd.subdivi_edge_no], gd);

	}
	
	//df_set_polyline_halfedge();
}


int curvenetwork::Patch::df_set_subpatch_over_hexagon(patchgen::SubPatchData *s, patchgen::GeometryData gd)
{

	if (s->num_sides == 5) {
		df_set_subpatch_pentagon(s, gd);
		return 0;
	}

	if (s->c_pointer_r->c_pointer_r != NULL) df_set_subpatch_over_hexagon(s->c_pointer_r, gd);
	if (s->c_pointer_l->c_pointer_r != NULL) df_set_subpatch_over_hexagon(s->c_pointer_l, gd);

	if (s->c_pointer_r->num_sides == 5) df_set_subpatch_pentagon(s->c_pointer_r, gd); 
	if (s->c_pointer_l->num_sides == 5) df_set_subpatch_pentagon(s->c_pointer_l, gd);

	std::vector<patchgen::SubPatchEdgeData *> share_edge;
	std::vector<std::vector<patchgen::SubPatchEdgeData *>> subdivi_edge;
	std::vector<std::vector<curvenetwork::Patch::VHandle>> subdivi_edge_vh;
	subdivi_edge_vh.resize(2);
	subdivi_edge_vh[0].resize(2);
	subdivi_edge_vh[1].resize(2);

	for (int i = 0; i < s->side_subdivi.size(); i++) {
		if (s->side_subdivi[i]->c_pointer_r != NULL && s->side_subdivi[i]->c_pointer_l != NULL) {
			std::vector<patchgen::SubPatchEdgeData *> tmp;
			tmp.push_back(s->side_subdivi[i]->c_pointer_r);
			tmp.push_back(s->side_subdivi[i]->c_pointer_l);
			subdivi_edge.push_back(tmp);
		}
	}
	for (int i = 0; i < s->c_pointer_r->side_subdivi.size(); i++) {
		if (s->c_pointer_r->side_subdivi[i]->p_pointer == NULL 
			&& s->c_pointer_r->side_subdivi[i]->same_edge_pointer == NULL) {

			subdivi_edge_vh[0][0] = (subdivi_edge[0][0]->v_handle_from.idx() == s->c_pointer_r->side_subdivi[i]->v_handle_to.idx()) 
				? subdivi_edge[0][0]->v_handle_from : subdivi_edge[0][0]->v_handle_to;
			subdivi_edge_vh[1][0] = (subdivi_edge[1][0]->v_handle_from.idx() == s->c_pointer_r->side_subdivi[i]->v_handle_to.idx())
				? subdivi_edge[1][0]->v_handle_from : subdivi_edge[1][0]->v_handle_to;
		}
	}
	for (int i = 0; i < s->c_pointer_l->side_subdivi.size(); i++) {
		if (s->c_pointer_l->side_subdivi[i]->p_pointer == NULL
			&& s->c_pointer_l->side_subdivi[i]->same_edge_pointer == NULL) {

			subdivi_edge_vh[0][1] = (subdivi_edge[0][1]->v_handle_from.idx() == s->c_pointer_l->side_subdivi[i]->v_handle_to.idx())
				? subdivi_edge[0][1]->v_handle_from : subdivi_edge[0][1]->v_handle_to;
			subdivi_edge_vh[1][1] = (subdivi_edge[1][1]->v_handle_from.idx() == s->c_pointer_l->side_subdivi[i]->v_handle_to.idx())
				? subdivi_edge[1][1]->v_handle_from : subdivi_edge[1][1]->v_handle_to;
		}
	}
	auto h_now_r = opposite_halfedge_handle(halfedge_handle(subdivi_edge_vh[0][0]));
	auto h_now_l = next_halfedge_handle(opposite_halfedge_handle(halfedge_handle(subdivi_edge_vh[0][1])));
	auto h_now_r_c = prev_halfedge_handle(h_now_r);
	curvenetwork::Patch::HHandle f_h;
	int r_l_flag = 0;
	for (;;) {
		auto v_r = to_vertex_handle(h_now_r_c);
		if (data(v_r).patchgen.corner_index != -1) {
			if (v_r.idx() != subdivi_edge_vh[1][0].idx()) r_l_flag = 1;
			break;
		}
		h_now_r_c = prev_halfedge_handle(h_now_r_c);
	}

	if (r_l_flag == 0){
		for (;;){
			auto v_r = to_vertex_handle(h_now_r);
			auto v_l = from_vertex_handle(h_now_l);
			auto prev_h_now_l = prev_halfedge_handle(h_now_l);
			auto next_h_now_r = next_halfedge_handle(h_now_r);

			auto nh_0 = new_edge(v_r, v_l);
			auto op_nh_0 = opposite_halfedge_handle(nh_0);

			set_next_halfedge_handle(h_now_r, nh_0);
			set_prev_halfedge_handle(nh_0, h_now_r);
			set_next_halfedge_handle(nh_0, h_now_l);
			set_prev_halfedge_handle(h_now_l, nh_0);

			set_next_halfedge_handle(op_nh_0, next_h_now_r);
			set_prev_halfedge_handle(next_h_now_r, op_nh_0);
			set_next_halfedge_handle(prev_h_now_l, op_nh_0);
			set_prev_halfedge_handle(op_nh_0, prev_h_now_l);

			collapse(nh_0);//rが消える？
			int v_r_idx = v_r.idx();
			//vertexとhalfedgeとfaceをdeleteする？
			f_h = h_now_r;
			if (v_r_idx == subdivi_edge_vh[1][0].idx()) break;
			h_now_r = prev_halfedge_handle(h_now_r);
			h_now_l = next_halfedge_handle(h_now_l);
		}
	}
	else if (r_l_flag == 1){
		h_now_r = next_halfedge_handle(h_now_r);
		h_now_l = prev_halfedge_handle(h_now_l);
		for (;;){
			auto v_r = from_vertex_handle(h_now_r);
			auto v_l = to_vertex_handle(h_now_l);

			auto next_h_now_l = next_halfedge_handle(h_now_l);
			auto prev_h_now_r = prev_halfedge_handle(h_now_r);

			auto nh_0 = new_edge(v_r, v_l);
			auto op_nh_0 = opposite_halfedge_handle(nh_0);

			set_prev_halfedge_handle(h_now_r, nh_0);
			set_next_halfedge_handle(nh_0, h_now_r);
			set_prev_halfedge_handle(nh_0, h_now_l);
			set_next_halfedge_handle(h_now_l, nh_0);

			set_prev_halfedge_handle(op_nh_0, prev_h_now_r);
			set_next_halfedge_handle(prev_h_now_r, op_nh_0);
			set_prev_halfedge_handle(next_h_now_l, op_nh_0);
			set_next_halfedge_handle(op_nh_0, next_h_now_l);

			collapse(nh_0);//rが消える？
			int v_r_idx = v_r.idx();
			data(v_r).patchgen.corner_index = -1;
			data(v_l).patchgen.corner_index = -1;
			//vertexとhalfedgeとfaceをdeleteする？
			f_h = h_now_r;
			if (v_r_idx == subdivi_edge_vh[1][0].idx()) break;
			h_now_r = next_halfedge_handle(h_now_r);
			h_now_l = prev_halfedge_handle(h_now_l);
		}
	}
	for (int i = 0; i < s->c_pointer_r->side_subdivi.size(); i++) {
		if (s->c_pointer_r->side_subdivi[i]->same_edge_pointer != NULL) {
			s->c_pointer_r->side_subdivi[i]->same_edge_pointer->v_handle_from = s->c_pointer_r->side_subdivi[i]->v_handle_from;
			s->c_pointer_r->side_subdivi[i]->same_edge_pointer->v_handle_to = s->c_pointer_r->side_subdivi[i]->v_handle_to;
		}
	}
	for (int i = 0; i < s->c_pointer_l->side_subdivi.size(); i++) {
		if (s->c_pointer_l->side_subdivi[i]->same_edge_pointer != NULL) {
			s->c_pointer_l->side_subdivi[i]->same_edge_pointer->v_handle_from = s->c_pointer_l->side_subdivi[i]->v_handle_from;
			s->c_pointer_l->side_subdivi[i]->same_edge_pointer->v_handle_to = s->c_pointer_l->side_subdivi[i]->v_handle_to;
		}
	}
	
	for (auto f : faces()){
		if (halfedge_handle(f).idx() != -1) {
			set_halfedge_handle(f, f_h);
		
		}
	}

	for (int i = 0; i < s->side_subdivi.size(); i++) {

		if (s->side_subdivi[i]->c_pointer_r != NULL && s->side_subdivi[i]->c_pointer_l != NULL) {

			curvenetwork::Patch::VHandle t;
			for (int j = 0; j < subdivi_edge_vh.size(); j++) {
				if (subdivi_edge_vh[j][1].idx() == s->side_subdivi[i]->c_pointer_r->v_handle_from.idx()) t = s->side_subdivi[i]->c_pointer_r->v_handle_from;
				else if(subdivi_edge_vh[j][1].idx() == s->side_subdivi[i]->c_pointer_r->v_handle_to.idx()) t = s->side_subdivi[i]->c_pointer_r->v_handle_to;
				else if (subdivi_edge_vh[j][1].idx() == s->side_subdivi[i]->c_pointer_l->v_handle_from.idx()) t = s->side_subdivi[i]->c_pointer_l->v_handle_from;
				else if (subdivi_edge_vh[j][1].idx() == s->side_subdivi[i]->c_pointer_l->v_handle_to.idx()) t = s->side_subdivi[i]->c_pointer_l->v_handle_to;
			}

			curvenetwork::Patch::HHandle bt;
			for (auto vh = cvih_iter(t); vh.is_valid(); ++vh) {
				//std::cout << "vh(HHandle) = " << vh.handle().idx() << "\n";
				for (auto f : faces()) {
					for (auto vv = cfh_iter(f); vv.is_valid(); ++vv) {
						if (vv.handle().idx() == vh.handle().idx()) {
							bt = vh.handle();
							break;
						}
						else if (vv.handle().idx() == opposite_halfedge_handle(vh.handle()).idx()) {
							bt = opposite_halfedge_handle(vh.handle());
							break;
						}
					}
				}
			}
			/*
			std::cout << "t(VHandle) = " << t.idx() << "\n";
			std::cout << "bt(HHandle) = " << bt.idx() << "\n";

			if (bt.idx() == -1){
				//df_debug_vertices_data();
				//df_debug_halfedges_data();
				//df_debug_faces_data();
				for (;;);
			}
			*/
			data(t).patchgen.corner_index = -1;
			curvenetwork::Patch::HHandle n;
			n = next_halfedge_handle(bt);
			for (;;) {
				if (data(from_vertex_handle(n)).patchgen.corner_index != -1) break;
				n = next_halfedge_handle(n);
			}

			curvenetwork::Patch::HHandle p;
			p = prev_halfedge_handle(bt);
			for (;;) {
				if (data(from_vertex_handle(p)).patchgen.corner_index != -1) break;
				p = prev_halfedge_handle(p);
			}

			s->side_subdivi[i]->v_handle_to = from_vertex_handle(n);
			s->side_subdivi[i]->v_handle_from = from_vertex_handle(p);

		}
	}

	return 0;
}

int curvenetwork::Patch::df_set_subpatch_pentagon(patchgen::SubPatchData *s, patchgen::GeometryData gd)
{
	int num_v_base = s->num_sides;

	std::vector<curvenetwork::Patch::VHandle> v_base(num_v_base);
	for (int vi = 0; vi < num_v_base; vi++)
	{
		v_base[vi] = add_vertex(curvenetwork::Patch::Point());
		data(v_base[vi]).patchgen.corner_index = vi;
		data(v_base[vi]).patchgen.irregular_index = -1;
	}

	std::vector<curvenetwork::Patch::VHandle>  face_v_base;
	for (int fvi = 0; fvi < num_v_base; fvi++)
		face_v_base.push_back(v_base[fvi]);
	s->sub_face = add_face(face_v_base);

	auto iv = new_vertex();
	data(iv).patchgen.irregular_index = 0;
	data(iv).patchgen.corner_index = -1;

	int count = 0;
	curvenetwork::Patch::HHandle iv_h_1st;
	for (auto h = cfh_iter(s->sub_face); h.is_valid(); ++h) {
		if (data(to_vertex_handle(h)).patchgen.corner_index != -1) {

			s->side_subdivi[count]->v_handle_from = from_vertex_handle(h);
			s->side_subdivi[count]->v_handle_to = to_vertex_handle(h);

			auto sv = new_vertex();
			data(sv).patchgen.irregular_index = -1;
			data(sv).patchgen.corner_index = -1;
			split_edge(edge_handle(h), sv);
			curvenetwork::Patch::HHandle nh = new_edge(iv, sv);

			if (count == 0)	iv_h_1st = nh;

			auto op_nh = opposite_halfedge_handle(nh);
			auto hh = h.handle();
			auto prev_hh = prev_halfedge_handle(hh);
			auto op_prev_hh = opposite_halfedge_handle(prev_hh);
			auto op_hh = opposite_halfedge_handle(hh);

			set_next_halfedge_handle(nh, op_prev_hh);
			set_prev_halfedge_handle(op_prev_hh, nh);
			set_prev_halfedge_handle(op_nh, op_hh);
			set_next_halfedge_handle(op_hh, op_nh);

			if (count == 0){
				set_prev_halfedge_handle(nh, op_nh);
				set_next_halfedge_handle(op_nh, nh);
			} else {
				auto next_next_op_prev_hh = next_halfedge_handle(next_halfedge_handle(op_prev_hh));
				set_prev_halfedge_handle(nh, next_next_op_prev_hh);
				set_next_halfedge_handle(next_next_op_prev_hh, nh);
				set_next_halfedge_handle(op_nh, iv_h_1st);
				set_prev_halfedge_handle(iv_h_1st, op_nh);
			}
			count++;
		}
		set_halfedge_handle(iv, iv_h_1st);
	}
	df_insert_polychord_pentagon(s->sub_face, gd, s);

	return 0;
}


void curvenetwork::Patch::df_insert_polychord_pentagon(curvenetwork::Patch::FHandle pf, patchgen::GeometryData gd, patchgen::SubPatchData *s)
{
	Eigen::VectorXi ip;
	ip = patchgen::df_calc_ilp_pentagon_polychord(s, gd);

	for (int i = 0; i < s->side_subdivi.size(); i++){
		s->side_subdivi[i]->polyline = opposite_halfedge_handle(halfedge_handle(s->side_subdivi[i]->v_handle_to));
		s->side_subdivi[i]->num_polychord = ip(i) - 1;

		df_insert_polychord(s->side_subdivi[i]->num_polychord, s->side_subdivi[i]->polyline);
	}
}

void curvenetwork::Patch::df_set_polyline_halfedge()
{
	polyline.clear();

	for (auto f : faces()) {
		for (auto fh = cfh_iter(f); fh.is_valid(); ++fh) {
			if (!data(fh.handle()).patchgen.polyline_lock){
				std::vector<curvenetwork::Patch::HHandle> h_polyline;
				h_polyline = df_find_polyline_halfedge(fh);
				polyline.push_back(h_polyline[0]);
				data(h_polyline[1]).patchgen.polyline_lock = true;
			}
		}
	}
}


void curvenetwork::Patch::df_insert_polychord(int insert_num, curvenetwork::Patch::HHandle ph)
{
	if (insert_num == 0) return;

	for (int c = 0; c < insert_num; c++) {

		auto h_start = opposite_halfedge_handle(ph);
		std::vector<curvenetwork::Patch::VHandle> iv;

		for (auto h = h_start;;) {
			auto h_next = opposite_halfedge_handle(next_halfedge_handle(next_halfedge_handle(h)));

			auto pf = point(from_vertex_handle(h));
			auto pt = point(to_vertex_handle(h));
			auto v_new = add_vertex((1 - 0.5) * pf + 0.5 * pt);
			iv.push_back(v_new);

			split_edge(edge_handle(h), v_new);

			if (!is_boundary(h)) break;
			if (h_next == h_start) break;

			h = h_next;
		}

		std::vector<curvenetwork::Patch::HHandle> ih;
		int ie_num = iv.size() - 1;
		for (int i = 0; i < ie_num; i++)
		{
			auto nh = new_edge(iv[i], iv[i + 1]);
			ih.push_back(nh);

			auto prev_hh_from_h = prev_halfedge_handle(halfedge_handle(from_vertex_handle(nh)));
			if (i > 0 && (prev_hh_from_h.idx() == ih[i - 1].idx()))
				prev_hh_from_h = opposite_halfedge_handle(halfedge_handle(from_vertex_handle(nh)));
			auto prev_prev_prev_hh_from_h = prev_halfedge_handle(prev_halfedge_handle(prev_hh_from_h));
			auto prev_prev_prev_prev_hh_from_h = prev_halfedge_handle(prev_halfedge_handle(prev_halfedge_handle(prev_hh_from_h)));
			auto next_prev_hh_from_h = next_halfedge_handle(prev_hh_from_h);
			auto hh_to_h = halfedge_handle(to_vertex_handle(nh));
			auto prev_hh_to_h = prev_halfedge_handle(hh_to_h);

			set_prev_halfedge_handle(nh, prev_hh_from_h);
			set_next_halfedge_handle(nh, prev_prev_prev_hh_from_h);
			set_prev_halfedge_handle(prev_prev_prev_hh_from_h, nh);
			set_next_halfedge_handle(prev_hh_from_h, nh);

			set_next_halfedge_handle(opposite_halfedge_handle(nh), next_prev_hh_from_h);
			set_prev_halfedge_handle(next_prev_hh_from_h, opposite_halfedge_handle(nh));
			set_prev_halfedge_handle(opposite_halfedge_handle(nh), prev_prev_prev_prev_hh_from_h);
			set_next_halfedge_handle(prev_prev_prev_prev_hh_from_h, opposite_halfedge_handle(nh));
		}
	}
}

void curvenetwork::Patch::df_draw_boundary(patchgen::GeometryData gd)
{
	auto draw_boundary = [&](GLenum mode) {
		int size = gd.num_sides;
		for (int i = 0; i < size; ++i) {
			glBegin(mode);
			for (int j = 0; j <= gd.side_subdivi[i]; ++j) {
				double t = i + j / static_cast<double>(gd.side_subdivi[i]);
				glVertex2d(patchgen::df_get_boundary_geometry(size, t, gd.corner_position));
			}
			glEnd();
		}
	};
	glColor3d(0, 0, 0);
	glLineWidth(4);
	draw_boundary(GL_LINE_STRIP);
	glPointSize(8);
	draw_boundary(GL_POINTS);
}

void curvenetwork::Patch::df_draw_selected_boundary(patchgen::GeometryData gd, int selected_side)
{
	glColor3d(1, 0, 0);
	auto draw_selected_boundary = [&](GLenum mode) {
		int size = gd.num_sides;
		glBegin(mode);
		for (int j = 0; j <= gd.side_subdivi[selected_side]; ++j) {
			double t = selected_side + j / static_cast<double>(gd.side_subdivi[selected_side]);
			glVertex2d(patchgen::df_get_boundary_geometry(size, t, gd.corner_position));
		}
		glEnd();
	};
	draw_selected_boundary(GL_LINE_STRIP);
	draw_selected_boundary(GL_POINTS);
}

void curvenetwork::Patch::df_draw_edge_subdivi(patchgen::GeometryData gd, Eigen::VectorXi change_subdivi)
{
	glColor3d(0, 0, 0);
	glPointSize(3);
	glLineWidth(3);
	int size = gd.num_sides;
	for (int i = 0; i < size; ++i) {
		Vector2d p0 = gd.corner_position[i];
		Vector2d p1 = gd.corner_position[(i + 1 >= size) ? 0 : i + 1];
		Vector2d tmp = p0 - p1;
		Vector2d t0;

		t0.x() = tmp.normalized().y() * 0.2;
		t0.y() = -tmp.normalized().x() * 0.2;
		Vector2d t1;
		t1.x() = -tmp.normalized().y() * 0.5;
		t1.y() = tmp.normalized().x() * 0.5;
		Vector2d q = 0.5 * (p0 + p1);
		glPushMatrix();
		glTranslated(q.x() + t0.x(), q.y() + t0.y(), 0);
		glScaled(0.001, 0.001, 1);

		string str = boost::lexical_cast<string>(change_subdivi[i]);

		if (change_subdivi[i] == gd.side_subdivi[i]) glColor3d(0, 0, 0);
		else glColor3d(1.0, 0, 0);

		for (char c : str)
			glutStrokeCharacter(GLUT_STROKE_MONO_ROMAN, c);
		glPopMatrix();
	}
}


void curvenetwork::Patch::df_set_subpatchdata(patchgen::GeometryData gd)
{

	
}
