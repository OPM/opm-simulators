#include <pybind11/pybind11.h>
#include <opm/simulators/flow/python/Pybind11Exporter.hpp>

void Opm::Pybind::export_all(py::module& m) {
    export_PyBlackOilSimulator(m);
}

PYBIND11_MODULE(simulators, m)
{
    Opm::Pybind::export_all(m);
}
