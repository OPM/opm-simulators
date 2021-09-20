#ifndef OPM_PYBIND11_EXPORTER_HEADER_INCLUDED
#define OPM_PYBIND11_EXPORTER_HEADER_INCLUDED

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
//#include <pybind11/embed.h>

namespace py = pybind11;

namespace Opm::Pybind {
    void export_all(py::module& m);
    void export_PyBlackOilSimulator(py::module& m);
}

#endif //OPM_PYBIND11_EXPORTER_HEADER_INCLUDED
