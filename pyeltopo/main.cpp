#include <memory>
#include <algorithm>
#include <iterator>
#include <mtao/types.h>
#include <fstream>
#include "eltopo.h"
#include <iostream>

#include <eltopo3d/eltopo.h>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;


PYBIND11_MODULE(pyeltopo, m) {
    py::class_<ElTopoTracker>(m, "ElTopoTracker")
        .def(py::init<const ElTopoTracker::CRefCV3d&, const ElTopoTracker::CRefCV3i&>())
        .def("get_triangles",&ElTopoTracker::get_triangles)
        .def("get_vertices",&ElTopoTracker::get_vertices)
        .def("integrate",&ElTopoTracker::integrate_py)
        .def("improve",&ElTopoTracker::improve)
        .def("step",&ElTopoTracker::step_py);
    m.def("make_tracker",[](py::EigenDRef<Eigen::MatrixXd>& a, py::EigenDRef<Eigen::MatrixXi>& b) {return make_tracker(a,b); });
}


