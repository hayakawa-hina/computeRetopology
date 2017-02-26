#pragma once
#include "ElemPtrT.hh"
#include <patchgen/PatchParam.h>
#include <patchgen/PatchVertexTraits.h>
#include <utility>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <kt84/openmesh/base/LaplaceDirect.hh>
#include <kt84/openmesh/base/LaplaceIterative.hh>
#include <kt84/openmesh/base/Utility.hh>
#ifndef NDEBUG
#   include <kt84/openmesh/base/DebugInfo.hh>
#endif
#include <kt84/graphics/DisplayList.hh>
#include <kt84/geometry/Polyline_PointNormal.hh>


#include <patchgen/GeometryData.h>
#include <patchgen/SubPatchData.h>
#include <patchgen/PatchHalfedgeTraits.h>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/System/config.h>
#include <OpenMesh/Core/Mesh/Status.hh>


namespace curvenetwork {
    struct PatchTraits : public OpenMesh::DefaultTraits {
        typedef OpenMesh::Vec3d Point;              // these are basically unused
        typedef OpenMesh::Vec3d Normal;
        
        VertexTraits
            , public patchgen::PatchVertexTraits
            , public kt84::LaplaceDirect_VertexTraits   <6>
            , public kt84::LaplaceIterative_VertexTraits<6>
        {
            kt84::PointNormal     pn;
            kt84::PointNormal     pn_temp;      // temporary data used by smoothing
            Eigen::Vector2d harmonic_uv;        // Canonical UVs used for SketchRetopo::compute_patch_interior_pn_harmonic()
            
            VertexT()
                : pn     (kt84::PointNormal::Zero())
                , pn_temp(kt84::PointNormal::Zero())
            {}
        };
        
        HalfedgeTraits
			, public patchgen::PatchHalfedgeTraits
            , public kt84::LaplaceDirect_HalfedgeTraits
            , public kt84::LaplaceIterative_HalfedgeTraits
        {
            Halfchain* halfchain;
            int        index_wrt_halfchain;
            
            HalfedgeT()
                : halfchain()
                , index_wrt_halfchain(-1)
            {}
        };
		//---------------------------------------------------
		VertexAttributes(OpenMesh::Attributes::Status);
		FaceAttributes(OpenMesh::Attributes::Status);
		EdgeAttributes(OpenMesh::Attributes::Status);
		//---------------------------------------------------
		
    };
    
    typedef OpenMesh::PolyMesh_ArrayKernelT<PatchTraits> PatchBase;
    
    struct Patch
        : public PatchBase
        , public kt84::LaplaceDirect   <PatchBase, Patch, 6>
        , public kt84::LaplaceIterative<PatchBase, Patch, 6>
        , public kt84::Utility         <PatchBase, Patch>
#ifndef NDEBUG
        , public kt84::DebugInfo       <PatchBase, Patch>
#endif
    {
        // basic data
        int  id;
        bool is_deleted;
        int  flag;
		explicit Patch(int id_ = -1);

		
        // connectivity info
        HalfchainPtr halfchain;
        
        // display list
        kt84::DisplayList displist_phong_fill;
        kt84::DisplayList displist_phong_line_interior;
        kt84::DisplayList displist_phong_line_boundary;
        kt84::DisplayList displist_flat_fill;
        kt84::DisplayList displist_flat_line_interior;
        kt84::DisplayList displist_flat_line_boundary;
        
        static double quadSize;                             // size of individual quads, used to determine Edgechain::num_subdiv. (don't know if it should really be declared here...)
        static bool   use_even_num_subdiv;
        static bool   prefer_rect_proc3;
        
        void clear();
        
        // connectivity-related functions
        void set_halfchain(Halfchain* halfchain_);
        HHandle opposite_boundary_halfedge(HHandle boundary_halfedge) const;
        Patch*  opposite_patch(HHandle boundary_halfedge) const;
        curvenetwork::Vertex* vertex_patch_to_curveNetwork(VHandle v) const;
        VHandle vertex_curveNetwork_to_patch(curvenetwork::Vertex* v) const;
        
        // other utility functions
        VHandle patch_corner         (int corner_index) const;
        std::pair<HHandle, HHandle> get_side_halfedges(int side_index) const;
        kt84::Polyline_PointNormal  get_side_polyline (int side_index) const;
        double distance(const Eigen::Vector3d& p) const;
        
        // static functions that only depends on Halfchain without requiring Patch instance (useful when instantiating Patch)
        static bool is_valid_for_patch(Halfchain*& halfchain);      // test if this halfchain forms a closed loop, and if so ensures that c->prev()->is_corner() holds
        static int num_corners(Halfchain* halfchain);
        static int num_subdiv(Halfchain* halfchain, int side_index);
        static Eigen::VectorXi num_subdiv(Halfchain* halfchain);
        static int     num_free_halfchains(Halfchain* halfchain, int side_index);
        static double  side_length        (Halfchain* halfchain, int side_index);
        static void    change_num_subdiv  (Halfchain* halfchain, int side_index, int target_num_subdiv);
        static void set_default_num_subdiv(Halfchain* halfchain);
        // non-static versions of above static functions with Patch::halfchain as arg
        int num_corners() const { return num_corners(halfchain); }
        int num_subdiv(int side_index) const { return num_subdiv(halfchain, side_index); }
        Eigen::VectorXi num_subdiv() const { return num_subdiv(&*halfchain); }
        
        // patch topology generation
        void add_padding(const patchgen::PatchParam& param);
       // patchgen::PatchParam generate_topology(const Eigen::VectorXi& num_subdiv);
        //void generate_topology(const patchgen::PatchParam& param);
        
        void set_halfedgeData ();
        void laplaceSolver_init();
        void set_boundary_pn();

        // rendering
        void render_phong_fill();
        void render_phong_line_interior();
        void render_phong_line_boundary();
        void render_flat_fill();
        void render_flat_line_interior();
        void render_flat_line_boundary();
        void render_irregular_vertices() const;
        void invalidate_displist();
		

		//debug
		void df_debug_vertices_data();
		void df_debug_faces_data();
		void df_debug_edges_data();
		void df_debug_halfedges_data();

		void df_polyline_seiretsu_kari5();

		std::vector<curvenetwork::Patch::HHandle> df_find_polyline_halfedge(curvenetwork::Patch::HHandle h);
		void df_set_geometry(patchgen::GeometryData gd);
		void df_set_base_patch(patchgen::GeometryData gd);
		void df_set_polyline_halfedge();
		void df_insert_polychord(int insert_num, curvenetwork::Patch::HHandle ph);
		void df_draw();
		void df_draw_boundary(patchgen::GeometryData gd);
		void df_draw_selected_boundary(patchgen::GeometryData gd, int selected_side);
		void df_draw_edge_subdivi(patchgen::GeometryData gd, Eigen::VectorXi change_subdivi);
		void df_set_subpatchdata(patchgen::GeometryData gd);
		int df_set_subpatch_over_hexagon(patchgen::SubPatchData *s, patchgen::GeometryData gd);
		int df_set_subpatch_pentagon(patchgen::SubPatchData *s, patchgen::GeometryData gd);
		void df_insert_polychord_pentagon(curvenetwork::Patch::FHandle pf, patchgen::GeometryData gd, patchgen::SubPatchData *s);
		//void draw_singularities();

		std::vector<HHandle> polyline;
	};
}
