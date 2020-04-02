#include "config.h"
#include <opm/parser/eclipse/Deck/Deck.hpp>
#include <opm/parser/eclipse/EclipseState/EclipseState.hpp>
#include <opm/parser/eclipse/EclipseState/Schedule/Schedule.hpp>
#include <opm/parser/eclipse/EclipseState/SummaryConfig/SummaryConfig.hpp>
#define OPM_FLOW_MAIN
#include <opm/simulators/flow/Main.hpp>
#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <pybind11/pybind11.h>
#include <pybind11/embed.h>
#include <iostream>

BEGIN_PROPERTIES

NEW_TYPE_TAG(EclFlowProblemMain);

END_PROPERTIES

namespace py = pybind11;

class BlackOilSimulator
{
public:
    BlackOilSimulator( const Opm::Deck& deck, 
                       const Opm::EclipseState& eclipseState, 
                       const Opm::Schedule& schedule,
                       const Opm::SummaryConfig& summaryConfig )
    {
        setDeck(deck);
        setEclipseState(eclipseState);
        setSchedule(schedule);
        setSummaryConfig(summaryConfig);
    }

    int run()
    {
        int argc = 1;
        char *argv[1] = {"flow"};
        using TypeTag = TTAG(EclFlowProblemMain);
        auto mainObject = Opm::Main<TypeTag>(
            argc, argv, deck_, eclipseState_, schedule_, summaryConfig_);
        return mainObject.run();
    }
    
    void setDeck( const Opm::Deck& deck )
    {
        deck_ = std::make_shared< Opm::Deck >(deck);
    }

    void setEclipseState( const Opm::EclipseState& eclipseState )
    {
        eclipseState_ = std::make_shared< Opm::EclipseState >(eclipseState);
    }

    void setSchedule( const Opm::Schedule& schedule )
    {
        schedule_ = std::make_shared< Opm::Schedule >(schedule);
    }

    void setSummaryConfig( const Opm::SummaryConfig& summaryConfig )
    {
        summaryConfig_ = std::make_shared< Opm::SummaryConfig >(summaryConfig);
    }

private:
    std::shared_ptr<Opm::Deck>          deck_;
    std::shared_ptr<Opm::EclipseState>  eclipseState_;
    std::shared_ptr<Opm::Schedule>      schedule_;
    std::shared_ptr<Opm::SummaryConfig> summaryConfig_;
};

PYBIND11_MODULE(simulators, m)
{
    py::class_<BlackOilSimulator>(m, "BlackOilSimulator")
        .def(py::init< const Opm::Deck&, const Opm::EclipseState&, const Opm::Schedule&, const Opm::SummaryConfig& >())
        .def("run", &BlackOilSimulator::run)
        .def("setDeck", &BlackOilSimulator::setDeck)
        .def("setEclipseState", &BlackOilSimulator::setEclipseState)
        .def("setSchedule", &BlackOilSimulator::setSchedule)
        .def("setSummaryConfig", &BlackOilSimulator::setSummaryConfig);
}
