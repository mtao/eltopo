#pragma once


//based off of ihttps://github.com/alecjacobson/gptoolbox/blob/master/mex/eltopo.cpp




#include <mtao/types.h>
#include <eltopo3d/surftrack.h>
#include <memory>

class SurfTrack;

class ElTopoTracker {
    public:
        using ColVectors3d = mtao::ColVectors<double,3>;
        using ColVectors3i = mtao::ColVectors<int,3>;
        using RefCV3d= Eigen::Ref<ColVectors3d>;
        using RefCV3i = Eigen::Ref<ColVectors3i>;
        using CRefCV3d= Eigen::Ref<const ColVectors3d>;
        using CRefCV3i = Eigen::Ref<const ColVectors3i>;
        ElTopoTracker& operator=(ElTopoTracker&&) = default;
        ~ElTopoTracker();
        ElTopoTracker(const CRefCV3d& V, const CRefCV3i& F, bool defrag_mesh = true, bool verbose = false);

        void split_edge(size_t edge_index);
        //split_triangle splits by the longest edge
        void split_triangle(size_t triangle_index);

        void initialize();
        void improve();
        double integrate(const CRefCV3d& Vnew, double dt);
        double step(const CRefCV3d& Vnew, double dt);
        double integrate_py(const Eigen::Ref<const Eigen::MatrixXd>& Vnew, double dt) {
            return integrate(Vnew,dt);
        }
        double step_py(const Eigen::Ref<const Eigen::MatrixXd>& Vnew, double dt) {
            return step(Vnew,dt);
        }

        void defrag_mesh();

        std::tuple<ColVectors3d,ColVectors3i> get_mesh() const;


        SurfTrackInitializationParameters& init_params() { return m_init_params; }
        ColVectors3d get_vertices() const;
        ColVectors3i get_triangles() const;
    private:
        std::unique_ptr<SurfTrack> m_surf = nullptr;
        SurfTrackInitializationParameters m_init_params;
        std::unique_ptr<SubdivisionScheme> m_subdivision_scheme;
        bool m_defrag_mesh = false;
        bool m_defrag_dirty = false;
        bool m_verbose = false;

};
