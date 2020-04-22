#include "config.h"
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#define FLOW_BLACKOIL_ONLY
#include <opm/simulators/flow/Main.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <iostream>
#include <string>

namespace py = pybind11;

class BlackOilSimulator
{
public:

    BlackOilSimulator( const std::string &deckFilename) : deckFilename_(deckFilename)
    {
    }

    int run()
    {
        auto mainObject = Opm::Main( deckFilename_ );
        return mainObject.runDynamic();
    }
    
private:
    const std::string deckFilename_;
};

PYBIND11_MODULE(simulators, m)
{
    py::class_<BlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init< const std::string& >())
        .def("run", &BlackOilSimulator::run);
}
