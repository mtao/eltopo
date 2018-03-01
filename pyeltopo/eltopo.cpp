#include <array>
#include <fstream>
#include <iterator>
#include <iostream>
#include <mtao/geometry/volume.h>
#include "eltopo.h"
#include <eltopo3d/subdivisionscheme.h>



ElTopoTracker::ElTopoTracker(const CRefCV3d& V, const CRefCV3i& F, bool defrag_mesh) {


    m_init_params.m_use_fraction = true;
    m_init_params.m_min_edge_length = .5;
    m_init_params.m_max_edge_length = 1.5;
    m_init_params.m_max_volume_change = .1;

    m_init_params.m_min_curvature_multiplier=1.0;
    m_init_params.m_max_curvature_multiplier=1.0;
    m_init_params.m_merge_proximity_epsilon=0.001;
    m_init_params.m_proximity_epsilon =1e-4;
    m_init_params.m_friction_coefficient=0.0;
    m_init_params.m_perform_improvement=true;
    m_init_params.m_allow_topology_changes=false;
    m_init_params.m_allow_non_manifold=false;
    m_init_params.m_collision_safety=true;


    m_subdivision_scheme.reset(new ButterflyScheme());
    m_init_params.m_subdivision_scheme=m_subdivision_scheme.get();


    std::vector<Vec3st> tris(F.cols());
    std::vector<Vec3d> verts(V.cols());
    std::vector<double> masses(V.cols());

    for(int i = 0; i < verts.size(); ++i) {
        Eigen::Map<mtao::Vector<double,3>> v(&verts[i][0]);
        v = V.col(i);
    }
    for(int i = 0; i < tris.size(); ++i) {
        Eigen::Map<mtao::Vector<size_t,3>> t(&tris[i][0]);
            t = F.col(i).cast<size_t>();
    }

    Eigen::Map<mtao::VectorX<double>>(masses.data(),V.cols()) = mtao::geometry::dual_volumes(V,F);


    m_surf = new SurfTrack(verts,tris,masses,m_init_params);


    if(defrag_mesh) {
        m_surf->defrag_mesh();
    }


    if ( m_surf->m_collision_safety )
    {
        m_surf->m_collision_pipeline.assert_mesh_is_intersection_free( false );      
    }

}

ElTopoTracker::~ElTopoTracker() {
    if(m_surf) {
        delete m_surf;
    }
}

auto ElTopoTracker::get_mesh() const -> std::tuple<ColVectors3d,ColVectors3i>{
    assert(m_surf);
    return std::make_tuple(get_vertices(),get_triangles());
}
auto ElTopoTracker::get_vertices() const ->ColVectors3d {
    assert(m_surf);
    auto&& pos = m_surf->get_positions();
    assert(pos.size() == m_surf->get_num_vertices());
    ColVectors3d V(3,pos.size());
    for(int i = 0; i < pos.size(); ++i) {
        V.col(i) = Eigen::Map<const mtao::Vector<double,3>>(&pos[i][0]);
    }
    return V;
}
auto ElTopoTracker::get_triangles() const -> ColVectors3i {
    assert(m_surf);
    auto&& tris = m_surf->m_mesh.get_triangles();
    ColVectors3i F(3,tris.size());
    for(int i = 0; i < tris.size(); ++i) {
        auto& t = tris[i];
        F.col(i) = Eigen::Map<const mtao::Vector<size_t,3>>(&t[0]).cast<int>();
    }
    assert(F.maxCoeff() < m_surf->get_num_vertices());
    return F;
}


void ElTopoTracker::improve() {
    assert(m_surf);
    // Improve
    m_surf->improve_mesh();
    // Topology changes
    m_surf->topology_changes();
    m_surf->defrag_mesh();
}


double ElTopoTracker::step(const CRefCV3d& V, double dt) {
    double r = integrate(V,dt);
    improve();
    return r;

}
double ElTopoTracker::integrate(const CRefCV3d& V, double dt) {
    assert(m_surf);

    double ret_dt;
    std::vector<Vec3d> verts(V.cols());
    for(int i = 0; i < verts.size(); ++i) {
        Eigen::Map<mtao::Vector<double,3>> v(&verts[i][0]);
        v = V.col(i);
    }

    m_surf->set_all_newpositions(verts);
    m_surf->integrate( dt, ret_dt );
    m_surf->set_positions_to_newpositions();
    return ret_dt;

}

